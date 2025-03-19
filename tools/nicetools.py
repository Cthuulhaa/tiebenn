import math
import os

import numpy as np
import pandas as pd

from obspy import UTCDateTime


def tt_theo_before_assoc(ev_time, teo_p_time, teo_s_time, pick, tol_p, tol_s):
    """

    This function decides if a predicted pick is close enough to the theoretical pick (for example, obtained using Taup)

    Args:
         ev_time (str) Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
         teo_p_time (float): Arrival time of the P-wave, in seconds after the event
         teo_s_time (float): Arrival time of the S-wave, in seconds after the event
         pick (SeisBench pick): The pick as predicted by the models in SeisBench
         tol_p (int or float): Seconds (in absolute value) in which a predicted P-pick is allowed to differ from the theoretical P-pick
         tol_s (int or float): Seconds (in absolute value) in which a predicted S-pick is allowed to differ from the theoretical S-pick
    """
    obs_time = pick.peak_time - ev_time

    if pick.phase == 'P':
       theo_obs_abs_diff = abs(obs_time - teo_p_time)
       if theo_obs_abs_diff > tol_p:
          verdict = 'fail'
       else:
            verdict = 'pass'
    else:
         theo_obs_abs_diff = abs(obs_time - teo_s_time)
         if theo_obs_abs_diff > tol_s:
            verdict = 'fail'
         else:
              verdict = 'pass'

    return verdict


def create_input_for_phassoc(outputs, seconds_before):
    """

    Slight modification of the outputs with phase picks, as necessary for the phase association

    Args:
         outputs (dict): The predicted phases/arrival times obtained from SeisBench for each station
         seconds_before (list): List of the seconds before the event
    Return:
           outputs_for_phassoc (dict): Modfied structure of the outputs dictionary
    """
    outputs_for_phassoc = {}

    for secs_bef in seconds_before:
        outputs_sec = {}
        for station in outputs:
            outputs_sec[station] = outputs[station][str(secs_bef)]
        outputs_for_phassoc[str(secs_bef)] = outputs_sec

    return outputs_for_phassoc


def generate_csv(streams, outputs, data, snr_data, ev_time):
    """

    Generates CSV files with pick detections

    Args:
         streams (obspy.stream.Stream): The retrieved ObsPy streams for each station
         outputs (dict): The predicted phases/arrival times obtained from SeisBench for each station
         data (dict): A dictionary with the nearest stations (ObsPy streams) to the event which contain useful channels for phase-picking
         snr_data (dict): Dictionary with the signal-to-noise ratio for each pick at each station
         ev_time (str): ev_time (str): Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
    Returns:
            prediction result files in CSV format for each station
    """
    ev_time = str(UTCDateTime(ev_time))
    os.mkdir(ev_time + '_tiebenn_loc/csv_picks')

    for stat in outputs[0]['station'].drop_duplicates():
        outputs_sta = outputs[0][outputs[0]['station'] == stat]

        net = []; sta = []; st_lat = []; st_lon = []; st_elv = []; phase = []; arrival_time = []; probability = []; snr = []
        for dets in range(len(outputs_sta)):
            net.append(data[stat]['network'])
            sta.append(stat)
            st_lat.append(data[stat]['coords'][0])
            st_lon.append(data[stat]['coords'][1])
            st_elv.append(data[stat]['coords'][2])

        for phs in outputs_sta['phase']:
            phase.append(phs)

            try:
                snr.append(snr_data[stat][phs])
            except Exception:
                   snr.append(np.nan)

        for arrt in outputs_sta['peak_time']:
            arrival_time.append(arrt)
        for prob in outputs_sta['peak_value']:
            probability.append(prob)

        df = {'network': net, 'station': sta, 'station_lat': st_lat, 'station_lon': st_lon, 'station_elv': st_elv, 'phase': phase, 'arrival_time': arrival_time, 'probability': probability, 'snr': snr}
        dataframe = pd.DataFrame(data=df)
        filename = ev_time + '_tiebenn_loc/csv_picks/' + stat + '_' + ev_time + '_prediction_results.csv'
        dataframe.to_csv(filename, index=False)

    return


def str2bool(v):
    """

    Tries to interpret an input as a boolean argument. Expected true boolean arguments are: yes, true, t, y, and 1. Expected false boolean arguments: no, false, f, n, 0

    Args:
         v (str): Argument to be interpreted as True or False
    Returns:
            v (bool): The string interpreted as boolean
    """
    if isinstance(v, bool):
       return v

    if v.lower() in ('yes', 'true', 't', 'y', '1'):
       return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
       return False
    else:
         raise argparse.ArgumentTypeError('Boolean value expected.')


def chan_comma(chanlist):
    """

    Transforms the channel list for a given station in a string of channels, each separated by a comma. This is the multiple channel format for a FDSN waveform request

    Args:
         chanlist (list): List of channels from which waveforms will be downloaded

    Returns:
         string (str): String of channels separated by a comma
    """
    string = ''

    for i in range(len(chanlist)):
        string = string + chanlist[i] + ','
    string = string.rstrip(string[-1])

    return string


def strmonth2num(strmonth):
    """

    It gets the month number corresponding to the 3-character month (e.g. for the input 'Jan' the output is 1)

    Args:
         strmonth (str): Month in 3-character format (not case-sensitive)

    Returns:
        nummonth (int): Month number
    """
    months = {'jan': 1, 'feb': 2, 'mar': 3, 'apr': 4, 'may': 5, 'jun': 6, 'jul': 7, 'aug': 8, 'sep': 9, 'oct': 10, 'nov': 11, 'dec': 12}

    nummonth = str(months[strmonth.lower()])

    return nummonth


def get_snr(data, picks, wlen):
    """

    Estimates the signal to noise ratio (SNR) of a seismogram. Adapted from the EQTransformer repository (https://github.com/smousavi05/EQTransformer)

    Args:
         data (ObsPy Stream): Stream with waveforms and metadata for a given station
         picks (dict): Sample point where a specific phase arrives
         wlen (int, positive): The length of the window in seconds for calculating the SNR

    Returns:
         snr (dict): A dictionary with the estimated SNR in db for each pick
    """
    snr = {}

    st_t = []; et_t = []; lendat = []
    for channels in range(len(data)):
        st_t.append(data[channels].stats.starttime)
        et_t.append(data[channels].stats.endtime)
        lendat.append(data[channels].stats.npts)

    starttime = min(st_t); endtime = max(et_t); lendata = max(lendat)

    def merge_pad(st, st_sta, st_end):
        cha = []; merge = False
        for ch in range(len(st)):
            if st[ch].stats.channel not in cha:
               cha.append(st[ch].stats.channel)
            else:
                merge = True

        if merge:
           st.merge()

        pad = False
        for ch in range(len(st)):
            if st[ch].stats.starttime - st_sta > 0 or st[ch].stats.endtime - st_end < 0:
               pad = True

        if pad:
           st.trim(starttime=st_sta, endtime=st_end, pad=True, fill_value=0)

        return st

    for p in picks['phase']:
        pick = picks[picks['phase'] == p]
        peak = UTCDateTime(str(pick['peak_time'].tolist()[0]))
        w_noise = data.copy()
        w_signal = data.copy()

        if (peak - starttime) >= wlen and (endtime - peak) > wlen:
           nw_sta = peak - wlen; nw_end = peak
           sw_sta = peak; sw_end = peak + wlen
        if (peak - starttime) < wlen and (endtime - peak) > wlen:
           nw_sta = starttime; nw_end = peak
           sw_sta = peak; sw_end = peak + wlen
        if (endtime - peak) <= wlen:
           nw_sta = peak - wlen; nw_end = peak
           sw_sta = peak; sw_end = endtime

        nw1 = w_noise.trim(starttime=nw_sta, endtime=nw_end)
        sw1 = w_signal.trim(starttime=sw_sta, endtime=sw_end)

        nw1 = merge_pad(st=nw1, st_sta=nw_sta, st_end=nw_end)
        sw1 = merge_pad(st=sw1, st_sta=sw_sta, st_end=sw_end)

        snr[pick['phase'].tolist()[0]] = round(10 * math.log10((np.percentile(sw1, 95) / np.percentile(nw1, 95)) ** 2), 1)

    return snr
