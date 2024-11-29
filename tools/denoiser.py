import obspy, os, joblib, pickle
from scipy.signal import stft, istft
from tensorflow.keras.models import load_model
import numpy as np

def load_obj(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


def preprocessing(data, dt=1.0, **kwargs):
    data = data - np.mean(data)

    trace = obspy.Trace(data=data, header=dict(delta=dt))

    try:
        if kwargs['decimation_factor']:
            trace.decimate(kwargs['decimation_factor'], no_filter=True)
            dt = trace.stats.delta
    except KeyError:
        pass

    try:
        trace.filter(**kwargs['filter'])
    except KeyError:
        pass

    try:
        trace.taper(**kwargs['taper'])
    except KeyError:
        pass

    return trace.data, dt


def merge_traces(stream: obspy.Stream, header: dict):
    """
    
    Merging traces in one obspy Stream. Note, altough obspy has a merging routine with good tests, this function is faster for the required method.
    In older versions "st_denoised.merge(method=1, fill_value="interpolate")" was used to merge traces of one stream, but this line was heavily time consuming.
    27.11.2024: I modified the original function because it used asyncio to merge traces in parallel. For some reason, this made the code crash inside Tiebenn.

    Args:
         stream (ObsPy Stream): contains all traces that will be merged
         header (dict): Dictionary that contains information for the stats of the merged stream. For more information about the stats have a look on the stats of an obspy trace.

    Returns:
         stream_out (ObsPy Stream): One obspy stream object instead of several single overlapping traces.
    """
    array_len = 0

    for trace in stream:
        array_len += trace.stats.npts

    data = np.zeros(array_len)
    start = 0

    for trace in stream:
        data[start:start + trace.stats.npts] = trace.data
        start += trace.stats.npts

    stream_out = obspy.Stream()
    trace_out = obspy.Trace(data=data, header=header)
    stream_out += trace_out

    return stream_out


def predict(model_filename, config_filename, data_list,  optimizer='adam'):
    """

    Function to predict data in data_list.

    Args:
         model_filename (str): Full filename path of model or checkpoints
         config_filename (str): Full filename path of config file
         data_list (list): List that contains numpy arrays for denoising
         optimizer (str): tensorflow optimizer for Model. Necessary if ckpt_model is True, default is adam

    Returns:
         recovered (array): Array that contains recoverd signal and noise. Has shape len(data_list) * ts_length * 2, where recovered[i, :, 0] denotes the recovered signal and recovered[i, :, 1] the recovered noise
         transform_list (array): Array that contains all relevant transformations. Has shape len(data_list) * shape_of_transformation * config['channels'] + 1. First channel contains transformation of original signal, second channel contains transformation of recovered signal and third channel of recovered noise
         
    """
    config = load_obj(config_filename)

    model_dae = load_model(model_filename)
    input_shape = (model_dae.input_shape[1], model_dae.input_shape[2])

    X = np.empty(shape=(len(data_list), *input_shape, config['channels']), dtype='float')
    transform_list = np.empty(shape=(len(data_list), *input_shape, config['channels'] + 1), dtype='complex')
    scales = []
    dj = []
    norm_factors = []
    dt = config['dt']
    mean_values = []

    if config['channels'] == 1:
        phases = []

    for i, array in enumerate(data_list):
        signal_tmp = array[:config['ts_length']]
        signal, dt = preprocessing(data=signal_tmp, dt=config['dt'], decimation_factor=config['decimation_factor'])
        mean_values.append(np.mean(signal))
        norm = np.max(np.abs(signal))
        if norm == 0:
            norm = 1
        signal = signal / norm
        norm_factors.append(norm)

        freqs, _, cns = stft(signal, fs=1 / dt, **config['kwargs'])

        if i == 0:
            recovered = np.empty(shape=(len(data_list), len(signal), 2), dtype='float')

        transform_list[i, :, :, 0] = cns

        X[i, :, :, 0] = cns.real / np.max(np.abs(cns.real))

        if config['channels'] == 1:
           phases.append(np.arctan2(cns.imag, cns.real))
        elif config['channels'] == 2:
             X[i, :, :, 1] = cns.imag / np.max(np.abs(cns.imag))
        else:
            msg = 'Channel number cannot exceed 2...'
            raise ValueError(msg)

    X_pred = model_dae.predict(X, verbose=0)

    for i in range(X_pred.shape[0]):
        if config['channels'] == 1:
           x_pred = X_pred[i, :, :, 0] * np.exp(1j * phases[i])
        elif config['channels'] == 2:
             transform_list[i, :, :, 1] = transform_list[i, :, :, 0] * X_pred[i, :, :, 0]
             transform_list[i, :, :, 2] = transform_list[i, :, :, 0] * X_pred[i, :, :, 1]

        t, rec_signal = istft(transform_list[i, :, :, 1], fs=1 / dt, **config['kwargs'])
        t, rec_noise = istft(transform_list[i, :, :, 2], fs=1 / dt, **config['kwargs'])

        rec_signal = np.real(rec_signal * norm_factors[i])
        rec_noise = np.real(rec_noise * norm_factors[i])

        rec_signal += mean_values[i]
        rec_noise += mean_values[i]

        recovered[i, :, 0] = rec_signal
        recovered[i, :, 1] = rec_noise

    return recovered, transform_list, freqs


def denoising_trace(trace, model_filename, config_filename, overlap=0.8, chunksize=None, verbose=True, **kwargs):
    """

    Denoising of an obspy Trace object using a trained Denoising Autoencoder of Tibi, as coded by J. Heuel (https://github.com/JanisHe/seisDAE)

    Args:
         trace (ObsPy Trace): obspy Trace to be denoised
         model_filename (str): full path to the trained denoising model
         config_filename (str): full path to the binary config file for the denoising model
         overlap (float): overlap between neighbouring elements in trace [0, 1]
         chunksize (int): for denosing of large traces, a trace is splitted into parts of chunksize, otherwise the data might not fit into memory.

    Returns:
         st_denoised, st_noise (ObsPy Trace): denoised and noisy traces, respectively
    """

    config = load_obj(config_filename)

    if trace.stats.delta != config['dt']:
       trace.resample(sampling_rate=1 / config['dt'])

    if config['ts_length'] > trace.stats.npts:
       len_conc = config['ts_length'] - trace.stats.npts
       conc = trace.data[-1] * np.ones(len_conc)
       trace.data = np.append(trace.data, conc)

    data_list = []
    starttime_list = []
    start = 0
    end = config['ts_length']

    while end <= trace.stats.npts:
          data_list.append(trace.data[start:end])
          starttime_list.append(trace.stats.starttime + start * trace.stats.delta)
          start += int(config['ts_length'] * (1 - overlap))
          end = start + config['ts_length']
    if end + 1 > trace.stats.npts:
        start = trace.stats.npts - config['ts_length']
        starttime_list.append(trace.stats.starttime + start * trace.stats.delta)
        data_list.append(trace.data[start:])

    if chunksize is not None and chunksize >= 1:
       chunks = int(np.ceil(len(data_list) / chunksize))

       for j in range(chunks):
           d = data_list[int(j * chunksize): int((j + 1) * chunksize)]
           r, _, _ = predict(model_filename, config_filename, d, **kwargs)

           if j == 0:
              recovered = copy.copy(r)
           else:
                recovered = np.concatenate((recovered, r))
    else:
         recovered, _, _ = predict(model_filename, config_filename, data_list, **kwargs)

    st_denoised = obspy.Stream()
    st_noise = obspy.Stream()

    if config['decimation_factor'] is not None:
        dt = trace.stats.delta * config['decimation_factor']
    else:
        dt = trace.stats.delta

    for i in range(recovered.shape[0]):
        tr_den = obspy.Trace(data=recovered[i, :, 0], header=dict(delta=dt, starttime=starttime_list[i], station=trace.stats.station, network=trace.stats.network, location=trace.stats.location, channel=trace.stats.channel))
        tr_noi = obspy.Trace(data=recovered[i, :, 1], header=dict(delta=dt, starttime=starttime_list[i], station=trace.stats.station, network=trace.stats.network, location=trace.stats.location, channel=trace.stats.channel))

        if overlap > 0:
            if i == 0:
                tr_den.trim(endtime=tr_den.stats.endtime - tr_den.stats.npts * tr_den.stats.delta * overlap / 2)
                tr_noi.trim(endtime=tr_noi.stats.endtime - tr_noi.stats.npts * tr_noi.stats.delta * overlap / 2)
            elif i == recovered.shape[0] - 1:
                tr_den.trim(starttime=st_denoised[-1].stats.endtime + dt)
                tr_noi.trim(starttime=st_noise[-1].stats.endtime + dt)
            else:
                tr_den.trim(starttime=st_denoised[-1].stats.endtime + dt, endtime=tr_den.stats.endtime - tr_den.stats.npts * tr_den.stats.delta * overlap / 2)
                tr_noi.trim(starttime=st_noise[-1].stats.endtime + dt, endtime=tr_noi.stats.endtime - tr_noi.stats.npts * tr_noi.stats.delta * overlap / 2)

        st_denoised += tr_den
        st_noise += tr_noi

    st_denoised = merge_traces(st_denoised, header=dict(delta=dt, starttime=starttime_list[0], station=trace.stats.station, network=trace.stats.network, location=trace.stats.location, channel=trace.stats.channel))

    st_noise = merge_traces(st_noise, header=dict(delta=dt, starttime=starttime_list[0], station=trace.stats.station, network=trace.stats.network, location=trace.stats.location, channel=trace.stats.channel))

    if verbose:
        print(f'Successfully denoised {trace.id} between {trace.stats.starttime} and {trace.stats.endtime}')

    return st_denoised[0], st_noise[0]


def denoising_stream(stream, model_filename, config_filename, overlap=0.8, chunksize=None, parallel=False, verbose=True, **kwargs):
    """

    Denoises an obspy stream and returns the recovered signal and noise as two separate streams.
    Note, the parameters not mentioned in the description are given in denoising_trace.

    :param stream: obspy stream object
    :param model_filename: Filename of a trained autoencoder
    :param config_filename: Filename of the config file that belongs to the autoencoder
    :param overlap: overlap between neighbouring segments. Default is 0.8
    :param chunksize: If the model has many parameters, chunksize splits all 60 s windows
                      into small chunks to reduce the memory. A value of 600 - 800 is recommended. Default is None.
    :param parallel: bool, default is False
                     If True, denoising is done in parallel otherwise one a single CPU
    :returns: stream of recovered signal and noise
    """
    if len(stream) == 0:
        msg = "The input stream does not contain any data.\n{}".format(str(stream))
        raise ValueError(msg)

    st_rec_signal = obspy.Stream()
    st_rec_noise = obspy.Stream()

    if parallel is True:
        if len(stream) <= int(os.cpu_count() / 2):
            n_jobs = len(stream)
        else:
            n_jobs = int(os.cpu_count() / 2)

        pool = joblib.Parallel(n_jobs=n_jobs, backend="multiprocessing", prefer="processes")
        out = pool(joblib.delayed(denoising_trace)(trace=trace, model_filename=model_filename, config_filename=config_filename, verbose=verbose, overlap=overlap, chunksize=chunksize, **kwargs) for trace in stream)

        for traces in out:
            st_rec_signal += traces[0]
            st_rec_noise += traces[1]
    else:
         for trace in stream:
             try:
                 tr_signal, tr_noise = denoising_trace(trace=trace, model_filename=model_filename, config_filename=config_filename, overlap=overlap, chunksize=chunksize, verbose=verbose, **kwargs)
                 st_rec_signal += tr_signal
                 st_rec_noise += tr_noise

             except:
                    pass

    return st_rec_signal, st_rec_noise
