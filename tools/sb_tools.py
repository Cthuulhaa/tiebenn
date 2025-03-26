import glob
import os

from joblib import Parallel, delayed
import numpy as np
import pandas as pd
import torch

import seisbench.models as sbm
from obspy import Stream, UTCDateTime
from obspy.clients.fdsn import Client as client_fdsn
from obspy.clients.filesystem.sds import Client as client_sds
from obspy.geodetics.base import gps2dist_azimuth

from .nicetools import chan_comma, create_input_for_phassoc, generate_csv, get_snr
from .visualization import plot_assoc
from .visualization import plotpicks_sb as plot
from .visualization import plotpicks_sb_mw as plot_mw


def picks_sb(ev_time, ev_lon, ev_lat, data, max_dist, client, picker, velmod, secs_before, phase_assoc, pick_sel, mult_windows, plotpicks, min_detections=3, denoise=True):
    """

    Produces probability functions and picks for P- and S-wave arrivals using SeisBench.

    Args:
         ev_time (str): Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
         ev_lon (float): Longitude of the located event
         ev_lat (float): Latitude of the located event
         data (dict): A dictionary with the nearest stations (ObsPy streams) to the event which contain useful channels for phase-picking
         max_dist (float): The stations used to seismic detection and P- and S-traveltime picks will be within this radius (in km)
         client (str): Choose between SDS or FDSN ObsPy client to retrieve waveforms. If SDS, make sure the SDS-directory is correctly defined
         picker (str): The pretrained model for phase picking: SeisBench_Phasenet or SeisBench_EQTransformer
         velmod (int): The velocity model used for phase association with PyOcto
         secs_before (float or list): If the picks will be predicted on multiple windows, this parameter is a list of seconds to be substracted from the event time. Otherwise, it is a constant value to be substracted from the event time
         phase_assoc (str): Implements a phase associator to the detected phase picks. Options are the GaMMA associator (https://github.com/AI4EPS/GaMMA) and PyOcto (https://github.com/yetinam/pyocto)
         pick_sel (char): If picks are calculated on multiple time windows, the selected criterion is applied for selecting picks. Options are: max_prob (for a given station, it looks for the pick(s) --P and/or S-- on that station at the time window where the probability was maximal), min_res (for a given station, it selects the pick(s) --P and/or S-- on time windows where the residual has the minimum value. This option is only available for the phase association with PyOcto)
         mult_windows (bool): If true, it will fetch waveforms for pick prediction with SeisBench starting at different times with respect to the event time
         plotpicks (bool): If True it will plot for each station with detections the event waveforms on the three components and their P- and S-arrivals. It will also plot the time series of the probabilities of an event, of a P- and a S-arrival. The phase association will also be plotted under an individual directory
         min_detections (int): For a given max_dist, it will define the minimum amount of stations on which the event was detected. If this minimum amount of detections is not achieved, the station radius will be increased and the detection process will be repeated. Value must be >= 3. Default = 3
         denoise (bool): If True it will apply the deep neural network DeepDenoiser (Zhu et al. 2019) as implemented on Seisbench (https://github.com/seisbench/seisbench) to denoise the retrieved waveforms before phase picking. Default is True

    Returns:
         <station>_<datetime>_prediction_results.csv: .csv files containing, for each station with sucessfully detected arrival times, all the detections and picking results.
         <station><datetime>.pdf (optional): .pdf files with plots of the waveforms, P- and S- probability functions and picks on each station
    """
    if phase_assoc.lower() == 'gamma' and pick_sel == 'min_res':
       pick_sel = 'max_prob'
       print('pick_sel set to max_prob')

    starttime = UTCDateTime(ev_time)

    if mult_windows == False:
       start_t = UTCDateTime(ev_time) - secs_before
       end_t = start_t + 60.

    data1 = {}
    data2 = {}
    data3 = {}

    for i1 in data:
        if float(data[i1]['epic_distance']) <= max_dist:
           data1[i1] = data[i1]
        if float(data[i1]['epic_distance']) <= max_dist + 100 and float(data[i1]['epic_distance']) > max_dist:
           data2[i1] = data[i1]
        if float(data[i1]['epic_distance']) <= max_dist + 200 and float(data[i1]['epic_distance']) > max_dist + 100:
           data3[i1] = data[i1]

    distances_dict = {}
    for et in data:
        distances_dict[et] = '{:.2f}'.format(gps2dist_azimuth(data[et]['coords'][0], data[et]['coords'][1], ev_lat, ev_lon)[0] * 0.001)

    detections = 0
    count = 1
    switch1 = 0; switch2 = 0; switch3 = 0

    while detections < min_detections:
         if count == 1:
            if len(data1) < min_detections:
               count = 2
               max_dist = max_dist + 100
            else:
                 if switch1 == 0:
                    print('ObsPy: Fetching station waveforms...')
                    if client.lower() in ('sds'):
                       sds_dir = '/SDS'
                       streams1 = get_streams_sds(start_time=starttime, data=data1, secs_before=secs_before, mult_windows=mult_windows, sds_dir=sds_dir)

                       to_delete = []
                       for i in data1.keys():
                           if i in streams1.keys():
                              to_delete.append(i)
                       data1_for_fdsn = data1.copy()

                       for i in to_delete:
                           data1_for_fdsn.pop(i)

                       print('Trying to include additional stations from FDSN server...')
                       streams1_fdsn = get_streams_fdsn_bulk(start_time=starttime, data=data1_for_fdsn, mult_windows=mult_windows, secs_before=secs_before)

                       if type(streams1_fdsn) != type(None):
                          streams1.update(streams1_fdsn)
                       else:
                            print('No additional streams retrieved...')

                    else:
                         streams1 = get_streams_fdsn_bulk(start_time=starttime, data=data1, mult_windows=mult_windows, secs_before=secs_before)
                    switch1 = 1

         if count == 2:
            if len(data1) + len(data2) < min_detections:
               count = 3
               max_dist = max_dist + 100
            else:
                 if switch1 == 0:
                    print('ObsPy: Fetching station waveforms...')
                    if client.lower() in ('sds'):
                       sds_dir = '/SDS'
                       streams1 = get_streams_sds(start_time=starttime, data=data1, secs_before=secs_before, mult_windows=mult_windows, sds_dir=sds_dir)

                       to_delete = []
                       for i in data1.keys():
                           if i in streams1.keys():
                              to_delete.append(i)
                       data1_for_fdsn = data1.copy()

                       for i in to_delete:
                           data1_for_fdsn.pop(i)

                       print('Trying to include additional stations from FDSN server...')
                       streams1_fdsn = get_streams_fdsn_bulk(start_time=starttime, data=data1_for_fdsn, mult_windows=mult_windows, secs_before=secs_before)

                       if type(streams1_fdsn) != type(None):
                          streams1.update(streams1_fdsn)
                       else:
                            print('No additional streams found...')
                    else:
                         streams1 = get_streams_fdsn_bulk(start_time=starttime, data=data1, mult_windows=mult_windows, secs_before=secs_before)
                    switch1 = 1
                 if switch2 == 0:
                    print('ObsPy: Fetching station waveforms...')
                    if client.lower() in ('sds'):
                       sds_dir = '/SDS'
                       streams2 = get_streams_sds(start_time=starttime, data=data2, secs_before=secs_before, mult_windows=mult_windows, sds_dir=sds_dir)

                       to_delete = []
                       for i in data2.keys():
                           if i in streams2.keys():
                              to_delete.append(i)
                       data2_for_fdsn = data2.copy()

                       for i in to_delete:
                           data2_for_fdsn.pop(i)

                       print('Trying to include additional stations from FDSN server...')
                       streams2_fdsn = get_streams_fdsn_bulk(start_time=starttime, data=data2_for_fdsn, mult_windows=mult_windows, secs_before=secs_before)

                       if type(streams2_fdsn) != type(None):
                          streams2.update(streams2_fdsn)
                       else:
                            print('No additional streams found...')
                    else:
                         streams2 = get_streams_fdsn_bulk(start_time=starttime, data=data2, mult_windows=mult_windows, secs_before=secs_before)
                    switch2 = 1

         if count == 3:
            if len(data1) + len(data2) + len(data3) < min_detections:
               print('STOP: Not enough stations with HH*, EH*, SH* or DN* channels available within ca 200 km...')
               break
            else:
                 if switch1 == 0:
                    print('ObsPy: Fetching station waveforms...')
                    if client.lower() in ('sds'):
                       sds_dir = '/SDS'

                       streams1 = get_streams_sds(start_time=starttime, data=data1, secs_before=secs_before, mult_windows=mult_windows, sds_dir=sds_dir)

                       to_delete = []
                       for i in data1.keys():
                           if i in streams1.keys():
                              to_delete.append(i)
                       data1_for_fdsn = data1.copy()

                       for i in to_delete:
                           data1_for_fdsn.pop(i)

                       print('Trying to include additional stations from FDSN server...')

                       streams1_fdsn = get_streams_fdsn_bulk(start_time=starttime, data=data1_for_fdsn, mult_windows=mult_windows, secs_before=secs_before)

                       if type(streams1_fdsn) != type(None):
                          streams1.update(streams1_fdsn)
                       else:
                            print('No additional streams found...')
                    else:
                        streams1 = get_streams_fdsn_bulk(start_time=starttime, data=data1, mult_windows=mult_windows, secs_before=secs_before)
                    switch1 = 1
                 if switch2 == 0:
                    print('ObsPy: Fetching station waveforms...')
                    if client.lower() in ('sds'):
                       sds_dir = '/SDS'
                       streams2 = get_streams_sds(start_time=starttime, data=data2, secs_before=secs_before, mult_windows=mult_windows, sds_dir=sds_dir)

                       to_delete = []
                       for i in data2.keys():
                           if i in streams2.keys():
                              to_delete.append(i)
                       data2_for_fdsn = data2.copy()

                       for i in to_delete:
                           data2_for_fdsn.pop(i)

                       print('Trying to include additional stations from FDSN server...')
                       streams2_fdsn = get_streams_fdsn_bulk(start_time=starttime, data=data2_for_fdsn, mult_windows=mult_windows, secs_before=secs_before)

                       if type(streams2_fdsn) != type(None):
                          streams2.update(streams2_fdsn)
                       else:
                            print('No additional streams found...')
                    else:
                        streams2 = get_streams_fdsn_bulk(start_time=starttime, data=data2, mult_windows=mult_windows, secs_before=secs_before)
                    switch2 = 1
                 if switch3 == 0:
                    print('ObsPy: Fetching station waveforms...')
                    if client.lower() in ('sds'):
                       sds_dir = '/SDS'
                       streams3 = get_streams_sds(start_time=starttime, data=data3, secs_before=secs_before, mult_windows=mult_windows, sds_dir=sds_dir)

                       to_delete = []
                       for i in data3.keys():
                           if i in streams3.keys():
                              to_delete.append(i)
                       data3_for_fdsn = data3.copy()

                       for i in to_delete:
                           data3_for_fdsn.pop(i)

                       print('Trying to include additional stations from FDSN server...')
                       streams3_fdsn = get_streams_fdsn_bulk(start_time=starttime, data=data3_for_fdsn, mult_windows=mult_windows, secs_before=secs_before)

                       if type(streams3_fdsn) != type(None):
                          streams3.update(streams3_fdsn)
                       else:
                            print('No additional streams found...')
                    else:
                        streams3 = get_streams_fdsn_bulk(start_time=starttime, data=data3, mult_windows=mult_windows, secs_before=secs_before)
                    switch3 = 1

         if 'streams1' in locals() and 'streams2' not in locals() and 'streams3' not in locals():
           if len(streams1) < min_detections:
              print('-------------------------------------------------------------------------------------------------------------------------------')
              print('WARNING: Less than', min_detections, 'stations available for P- and S-wave picking within', max_dist, 'km. Increasing radius...')
              print('-------------------------------------------------------------------------------------------------------------------------------')
              count = count + 1
              max_dist = max_dist + 100
              if max_dist > 200:
                 print('200 km epicentral distance limit reached. Stop.')
                 break
              else:
                   continue
           else:
                streams = streams1
         elif 'streams1' in locals() and 'streams2' in locals() and 'streams3' not in locals():
             if len(streams1) + len(streams2) < min_detections:
               print('-------------------------------------------------------------------------------------------------------------------------------')
               print('WARNING: Less than', min_detections, 'stations available for P- and S-wave picking within', max_dist, 'km. Increasing radius...')
               print('-------------------------------------------------------------------------------------------------------------------------------')
               if max_dist + 100 < 200:
                  max_dist = max_dist + 100
                  count = count + 1
                  continue
               elif max_dist < 200 and max_dist + 100 > 200:
                    max_dist = 200
                    count = count + 1
                    continue
               else:
                    print('200 km epicentral distance limit reached. Stop.')
                    break
             else:
                  streams = {**streams1, **streams2}
         elif 'streams1' in locals() and 'streams2' in locals() and 'streams3' in locals():
             if len(streams1) + len(streams2) + len(streams3) < min_detections:
                print('-------------------------------------------------------------------------------------------------------------------------------')
                print('WARNING: Less than', min_detections, 'stations available for P- and S-wave picking within', max_dist, 'km. Stop')
                print('-------------------------------------------------------------------------------------------------------------------------------')
                break
             else:
                  streams = {**streams1, **streams2, **streams3}

         print('**********************************************************************************************************')
         print('-#-#-#- Trying detections within', max_dist, 'km from epicenter... -#-#-#-')
         print('**********************************************************************************************************')

         if picker.lower() in ['sb_eqt', 'sb_eqtransformer', 'seisbench_eqt', 'seisbench_eqtransformer']:
            print('SeisBench: Creating EQTransformer picker from pre-trained model...')
            model = sbm.EQTransformer.from_pretrained('original') #XXX NOTE: Benchmark datasets which were used to train the EQT and PhaseNet models include: original, ethz, instance, scedc, stead, geofon, neic
         elif picker.lower() in ['sb_pn', 'sb_phasenet', 'seisbench_pn', 'seisbench_phasenet']:
              print('SeisBench: Creating PhaseNet picker from pre-trained model...')
              model = sbm.PhaseNet.from_pretrained('instance')

         if mult_windows:
            streams_for_plot = {}

         predictions = {}
         outputs = {}
         number_of_picks = 0

         n_jobs = calculate_njobs(streams)

         streams_preproc = Parallel(n_jobs=n_jobs, backend='threading')(delayed(waveform_preprocessing)(stream=streams[raw_stream], denoise=denoise, distance=distances_dict[raw_stream]) for raw_stream in streams)

         for s in streams_preproc:
             s_station = s[0].stats.station
             streams[s_station] = s

             if mult_windows:
                streams_for_plot[s_station] = s

         if not mult_windows:
            phasepicks = Parallel(n_jobs=n_jobs, backend='threading')(delayed(picks_singwin)(stream=streams[sta_], picker=picker, model=model) for sta_ in streams)

            for picks_ in phasepicks:
                s_station = picks_[0]
                predictions[s_station] = picks_[1]
                outputs[s_station] = picks_[2]

            for sta in outputs:
                if len(outputs[sta].picks) > 0:
                   number_of_picks = number_of_picks + 1

         else:
              phasepicks = parallel_phase_picking(station_streams=streams, starttime=starttime, picker=picker, model=model, start_list=secs_before)

              for picks_ in phasepicks:
                  s_station = picks_[0]
                  s_bef = str(picks_[1])

                  if s_station in predictions:
                     predictions[s_station][s_bef] = picks_[2]
                  else:
                       predictions[s_station] = {}
                       predictions[s_station][s_bef] = picks_[2]

                  if s_station in outputs:
                     outputs[s_station][s_bef] = picks_[3]
                  else:
                       outputs[s_station] = {}
                       outputs[s_station][s_bef] = picks_[3]

              for sta in outputs:
                  for secs in outputs[sta]:
                      if len(outputs[sta][secs].picks) != 0:
                         number_of_picks = number_of_picks + 1
                         break

         if number_of_picks >= min_detections:
            print('SeisBench: Automatic picks were obtained on %s stations' % str(number_of_picks))
         else:
              if max_dist >= 200:
                 print('STOP: The event was not detected on the minimum requested stations...')
                 break
                 return
              else:
                   print('P-/S-wave picks detected on less than', min_detections, 'stations within', max_dist, 'km. Increasing radius...')
                   count = count + 1
                   max_dist = max_dist + 100
                   continue

         max_num_stations = 80

         if len(outputs) > max_num_stations:
            sorted_distances = []

            for st in outputs:
                sorted_distances.append(float(data[st]['epic_distance']))
            sorted_distances.sort()
            sorted_distances = sorted_distances[0:max_num_stations]

            st_to_delete = []
            for st in outputs:
                if float(data[st]['epic_distance']) not in sorted_distances:
                   st_to_delete.append(st)

            for delete in st_to_delete:
                outputs.pop(delete)

         if phase_assoc.lower() == 'pyocto':
            from .pyoctools import phase_association
         else:
              from .gammaasoc import phase_association

         if mult_windows == False:
            try:
                print('Phase association of picks...')

                outputs_assoc = phase_association(outputs=outputs, data=data, velmod=velmod, ev_lon=ev_lon, ev_lat=ev_lat, ev_time=starttime, max_dist=max_dist, plot=plotpicks, mult_windows=mult_windows, secs_before=secs_before)
            except Exception as e:
                   print('Error in phase association:\n', e)
                   break

            if outputs_assoc == {}:
               print('STOP: not enough picks passed the comparison with theoretical arrival times! Skipping event.')
               break

            if len(outputs_assoc) == 1:
               outputs_stations_assoc = outputs_assoc[0]['station']

               number_of_picks = len(outputs_stations_assoc.drop_duplicates())
            else:
                 # XXX NOTE: In the following step, if the predicted picks were associated to more than one event, the event will be selected which has its event time closest to the original ev_time argument. The picks associated to this event will be used for depth estimation with NonLinLoc. An approach on how to proceed with the rest of events is yet to be decided. As more location/depth estimation examples are tested, it will be easier to figure out how to handle more complex cases, such as multiple events detected by PyOcto.
                 diff_evtime = 999999999
                 da_event = 999
                 outputs_assoc_one = {}
                 for i in outputs_assoc:
                     if abs(UTCDateTime(str(outputs_assoc[i]['time'].to_list()[0]).split('+')[0]) - starttime) < diff_evtime:
                        diff_evtime = abs(UTCDateTime(str(outputs_assoc[i]['time'].to_list()[0]).split('+')[0]) - starttime)
                        da_event = i
                 outputs_assoc_one[0] = outputs_assoc[da_event]
                 outputs_assoc = outputs_assoc_one

                 number_of_picks = len(outputs_assoc[0]['station'].drop_duplicates())
         else:
              outputs_for_phassoc = create_input_for_phassoc(outputs, seconds_before=secs_before)

              secs_bef_ = []
              for secs_bef in outputs_for_phassoc:
                  testing_lenght = 0
                  for test in outputs_for_phassoc[secs_bef]:
                      testing_lenght = testing_lenght + len(outputs_for_phassoc[secs_bef][test].picks)

                  if testing_lenght == 0:
                     print('No picks were produced for this event in the current time window. Skipping...')
                  else:
                       secs_bef_.append(secs_bef)

              print('Phase association of picks...')

              n_jobs = calculate_njobs(secs_bef_)

              associations = Parallel(n_jobs=n_jobs, backend='threading')(delayed(phase_association)(outputs=outputs_for_phassoc[s_bef], data=data, velmod=velmod, ev_lon=ev_lon, ev_lat=ev_lat, ev_time=starttime, plot=plotpicks, max_dist=max_dist,mult_windows=mult_windows, secs_before=s_bef) for s_bef in secs_bef_)

              outputs_assoc = {}
              for sbef in associations:
                  if sbef != None:
                     outputs_assoc[sbef[1]] = sbef[0]

              if len(outputs_assoc) == 0:
                 print('The predicted P- and S- picks for all the time windows were associated to no event. Skipping to next event in file...')
                 break

              for assocs in outputs_assoc:
                  if len(outputs_assoc[assocs]) > 1: # XXX NOTE: See comment above
                     diff_evtime = 999999999
                     da_event = 999
                     outputs_assoc_one = {}
                     for i in outputs_assoc[assocs]:
                         if abs(UTCDateTime(str(outputs_assoc[assocs][i]['time'].to_list()[0]).split('+')[0]) - starttime) < diff_evtime:
                            diff_evtime = abs(UTCDateTime(str(outputs_assoc[assocs][i]['time'].to_list()[0]).split('+')[0]) - starttime)
                            da_event = i
                     outputs_assoc_one[0] = outputs_assoc[assocs][da_event]
                     outputs_assoc[assocs] = outputs_assoc_one

              outputs_final = select_picks(outputs=outputs, outputs_assoc=outputs_assoc, phase_assoc=phase_assoc, mode=pick_sel)

              number_of_picks = len(outputs_final['station'].drop_duplicates())

         if number_of_picks >= min_detections:
            print('Picks on %s stations were associated to the event' % str(number_of_picks))
            detections = number_of_picks
         else:
              if max_dist >= 200:
                 print('STOP: Picks not associated on the minimum requested stations...')
                 break
                 return
              else:
                   print('P-/S-wave picks associated on less than', min_detections, 'stations. Increasing radius...')
                   count = count + 1
                   max_dist = max_dist + 100
                   continue

    if mult_windows == True:
       outputs_assoc = {}
       outputs_assoc[0] = outputs_final

    phases_snr = {}
    for sn in outputs_assoc[0]['station'].drop_duplicates().tolist():
        phases_snr[sn] = get_snr(data=streams[sn], picks=outputs_assoc[0][outputs_assoc[0]['station']==sn], wlen=2)

    os.mkdir(str(starttime) + '_tiebenn_loc')

    generate_csv(streams=streams, outputs=outputs_assoc, data=data, snr_data=phases_snr, ev_time=ev_time)

    if plotpicks:
       n_jobs = calculate_njobs(streams)

       try:
           os.mkdir(f'{starttime}_tiebenn_loc/plot_phase_association')
           for i in glob.glob('PhAssoc_*.pdf'):
               os.rename(i, f'{starttime}_tiebenn_loc/plot_phase_association/{i}')
       except:
              pass

       print('Generating waveform and plotpicks figures...')

       if not mult_windows:
          Parallel(n_jobs=n_jobs, backend='threading')(delayed(plot)(data=data, streams=streams[station], starttime=starttime, predictions=predictions, picks=outputs) for station in streams)
       else:
            Parallel(n_jobs=n_jobs, backend='threading')(delayed(plot_mw)(data=data, streams=streams_for_plot[station], starttime=starttime, predictions=predictions, picks=outputs, picks_final=outputs_final) for station in streams)

    if mult_windows == True:
       return streams_for_plot
    else:
         return streams


def deepdenoiser(stream):
    """

    DeepDenoiser (Zhu et al. 2019) as implemented on Seisbench.

    Args:
         stream (stream): Obspy stream with the waveforms to be denoised

    Returns:
            stream_denoise (stream): The stream with denoised waveforms
    """
    denoise = sbm.DeepDenoiser.from_pretrained('original')

    stream_denoised = denoise.annotate(stream)

    for ch in range(len(stream_denoised)):
        stream_denoised[ch].stats.channel = stream_denoised[ch].stats.channel.split('_')[1]

    return stream_denoised


def calculate_njobs(loop):
    """

    Calculates the number of parallel jobs to be run at parallelizations

    Args:
         loop: The variable which will determine the number of simultaneous jobs. It can be e.g. the number of stations
    """
    if len(loop) <= int(os.cpu_count() / 2):
       n_jobs = len(loop)
    else:
         n_jobs = int(os.cpu_count() / 2)

    return n_jobs


def process_client_waveforms(client, data, start_t, end_t):
    """
    Process waveform retrieval for a single client. It builds a bulk request for all stations associated with the client and retrieves the corresponding waveforms.

    Args:
         client (str): ObsPy client to request waveform bulk
         data (dict): A dictionary with information about the stations. Example: {'LANDS': {'network': 'SX', 'channels': ['HHN', 'HHE', 'HHZ'], 'coords': [51.526, 12.163, 115.0], 'client': 'BGR', 'epic_distance': '53.61'}
         start_t (str): Start of requested time window. Format: yyyy-mm-dd hh:mm:ss.ss
         end_t (str): End of requested time window. Format: yyyy-mm-dd hh:mm:ss.ss

    Returns:
            stream (dict): A dictionary with the retrieved streams which may be empty if no data is retrieved
    """
    bulk = []
    detect_duplicates = []

    for st in data:
        channels = chan_comma(data[st]['channels'])
        if data[st]['client'] == client:
            if data[st]['coords'] in detect_duplicates:
                print(f'Station duplicate detected for client {client}. Skipping station {st}...')
            else:
                bulk.append((data[st]['network'], st, '*', channels, start_t, end_t))
                detect_duplicates.append(data[st]['coords'])

    stream = Stream()
    if bulk:
        try:
            if client == 'BGR':
               try:
                   stream = client_fdsn('http://192.168.11.220:8080').get_waveforms_bulk(bulk)
               except:
                      stream = client_fdsn(client).get_waveforms_bulk(bulk)
            else:
                stream = client_fdsn(client).get_waveforms_bulk(bulk)
        except Exception as e:
            print(f'No waveforms retrieved from client {client}')

    return stream


def get_streams_fdsn_bulk(start_time, data, mult_windows, secs_before):
    """
    Uses ObsPy FDSN client to generate a bulk request to download streams with 60 second long waveforms for a set of stations.

    Args:
         start_time (UTCDateTime): Start of requested time window.
         data (dict): Station information dictionary.
         mult_windows (bool): If True, multiple time windows are used.
         secs_before (float or list): Time to subtract from start_time when mult_windows is False.

    Returns:
         stream_output (dict): Dictionary with the retrieved station streams.
    """
    if not mult_windows:
        start_t = start_time - secs_before
        end_t = start_t + 60
    else:
        start_t = start_time - 10
        end_t = start_time + 60

    client_list = ['BGR', 'LMU', 'GFZ', 'ODC', 'RASPISHAKE', 'RESIF', 'ETH', 'INGV', 'IPGP', 'NIEP', 'ORFEUS']

    n_jobs = calculate_njobs(client_list)

    streams_list = Parallel(n_jobs=n_jobs, backend='threading')(
        delayed(process_client_waveforms)(client, data, start_t, end_t)
        for client in client_list)

    combined_streams = Stream()
    for s in streams_list:
        combined_streams += s

    stream_output = {}
    if len(combined_streams) != 0:
        stations = []
        for tr in combined_streams:
            if tr.stats.station not in stations:
                stations.append(tr.stats.station)
        for station in stations:
            stream_out = combined_streams.select(station=station)
            if len(stream_out.get_gaps()) > 0:
                stream_out.merge()
            for tr in list(stream_out):
                if len(tr.data) == 0:
                    stream_out.remove(tr)
            if len(stream_out) > 0:
                try:
                    stream_output[station] = stream_out.trim(starttime=start_t, endtime=end_t)
                except Exception as e:
                    print(f'Error trimming stream for station {station}')
                    stream_output[station] = stream_out

    return stream_output


def get_streams_sds(start_time, data, secs_before, mult_windows, sds_dir):
    """

    Uses ObsPy SDS client to download streams with 60 second long waveforms for a set of stations.

    Args:
         start_time (str): Start of requested time window. Format: yyyy-mm-dd hh:mm:ss.ss
         data (dict): A dictionary with information about the stations. Example: {'LANDS': {'network': 'SX', 'channels': ['HHN', 'HHE', 'HHZ'], 'coords': [51.526, 12.163, 115.0], 'client': 'BGR', 'epic_distance': '53.61'}
         secs_before (float or list): If the picks will be predicted on multiple windows, this parameter is a list of seconds to be substracted from the event time. Otherwise, it is a constant value to be substracted from the event time
         mult_windows (bools): If True, the retrieved streams will be for multiple time windows around the event time
         sds_dir (str): Root directory of SDS archive

    Returns:
         streams (dict): A dictionary with the retrieved station streams
    """
    sds = client_sds(sds_dir)
    streams = {}

    if mult_windows == False:
       start_t = start_time - secs_before
       end_t = start_t + 60

    else:
         start_t = start_time - 10
         end_t = start_time + 60

    for st in data:
        print('Fetching data for station', st)
        stream = Stream()
        for ch in data[st]['channels']:
            trace = sds.get_waveforms(data[st]['network'], st, '*', channel=ch, starttime=start_t, endtime=end_t)

            if len(trace) != 0:
               if len(trace[0].data) > 0:
                  stream += trace

        if len(stream) != 0:
           if len(stream.get_gaps()) > 0:
              stream.merge()
           streams[st] = stream
        else:
             print('Fetching failed for station', st)

    return streams


def waveform_preprocessing(stream, denoise, distance):
    """

    Waveform processing with ObsPy

    Args:
         stream (ObsPy Stream): Stream with station channels to be pre-processed. Dead channels are removed, masked channels are split, data are detrended, bandpass-filtered and tapered
    Return:
         stream (ObsPy stream): Pre-processed ObsPy stream
    """

    station = stream[0].stats.station

    for if0 in stream:
        if len(if0.data) == 0: # We remove empty channels
           stream.remove(if0)

        masked = False # We look for masked data
        for mask in range(len(stream)):
            tr = stream[mask].data
            if np.ma.is_masked(tr) == True:
               masked = True

        if masked: # If masked channels, we split
           stream = stream.split()

    stream.detrend().filter('bandpass', freqmin=1, freqmax=25, corners=2, zerophase=True).taper(max_percentage=0.001, type='cosine', max_length=2)

    if denoise:
       crit_distance = 100.

       if float(distance) <= crit_distance:
          stream = deepdenoiser(stream)

    return stream


def picks_singwin(stream, picker, model):
    """

    Phase picks on stations for single-windows.

    Args:
         stream (ObsPy Stream): ObsPy Stream where the phases will be detected
         picker (str): Name of the phase picker to be used: either EQTransformer or PhaseNet
         model (SeisBench PhasePicker): The actual SeisBench model
    Returns:
         predictions (SeisBench Predictions): Probability time series of P-, S, and noise (for PhaseNet) or detection (for EQTransformer)
         outputs (SeisBench Outputs): Pick information. Time stamp, maximum probability, peak time, etc
    """
    station = stream[0].stats.station

    print('Phase picking for station %s' % station)

    if picker.lower() in ['sb_eqt', 'sb_eqtransformer', 'seisbench_eqt', 'seisbench_eqtransformer']:
       predictions = model.annotate(stream, overlap=3000, detection_threshold=0.25, P_threshold=0.2, S_threshold=0.15)
       outputs = model.classify(stream, overlap=3000, detection_threshold=0.25, P_threshold=0.2, S_threshold=0.15)

    elif picker.lower() in ['sb_pn', 'sb_phasenet', 'seisbench_pn', 'seisbench_phasenet']:
         predictions = model.annotate(stream, overlap=2800, P_threshold=0.2, S_threshold=0.15)
         outputs = model.classify(stream, overlap=2800, P_threshold=0.2, S_threshold=0.15)

    return station, predictions, outputs


def picks_mulwin(streams, start, starttime, picker, model):
    """

    Phase picks on stations for single-windows.

    Args:
         streams (ObsPy Stream): ObsPy Stream where the phases will be detected
         start (int or float): Seconds before event time for waveform time window to start
         starttime (ObsPy UTCDateTime): Event start time
         picker (str): Name of the phase picker to be used: either EQTransformer or PhaseNet
         model (SeisBench PhasePicker): The actual SeisBench model
    Returns:
         predictions (SeisBench Predictions): Probability time series of P-, S, and noise (for PhaseNet) or detection (for EQTransformer)
         outputs (SeisBench Outputs): Pick information. Time stamp, maximum probability, peak time, etc
    """
    station = streams[0].stats.station
    sta = starttime - start

    stream = Stream()
    dummy = streams.copy()

    for ch in range(len(streams)):
        stream += dummy[ch].trim(starttime=sta, endtime=sta+60)

    print('Phase picking for station', station, 'Starttime:', start, 'seconds before event time')

    if picker.lower() in ['sb_eqt', 'sb_eqtransformer', 'seisbench_eqt', 'seisbench_eqtransformer']:
       predictions = model.annotate(stream, overlap=3000, detection_threshold=0.25, P_threshold=0.2, S_threshold=0.15)
       outputs = model.classify(stream, overlap=3000, detection_threshold=0.25, P_threshold=0.2, S_threshold=0.15)

    elif picker.lower() in ['sb_pn', 'sb_phasenet', 'seisbench_pn', 'seisbench_phasenet']:
         predictions = model.annotate(stream, overlap=2800, P_threshold=0.2, S_threshold=0.15)
         outputs = model.classify(stream, overlap=2800, P_threshold=0.2, S_threshold=0.15)

    return predictions, outputs, start


def parallel_phase_picking(station_streams, starttime, picker, model, start_list):
    """

    Run phase picking in parallel for each station (stream) and each time window alignment.

    Args:
         station_streams (dict): Keys are station names, values are ObsPy Stream objects
         starttime (UTCDateTime): The event time
         picker (str): Name of the phase picker to use
         model: SeisBench phase picking model
         start_list (list): List of seconds before event time to try (e.g., [0, 2, 5, 10])

    Returns:
         results (list): A list of tuples (station, start, predictions, outputs)
    """
    n_jobs = calculate_njobs(np.zeros(len(station_streams) + len(start_list)))

    results = Parallel(n_jobs=n_jobs, backend='threading')(delayed(lambda st, start: (st[0], start, *picks_mulwin(st[1], start, starttime, picker, model)))(stream_item, start)
        for stream_item in station_streams.items()
        for start in start_list)

    return results


def select_picks(outputs, outputs_assoc, phase_assoc, mode):
    """

    When the multiple-window approach (Park et al. 2023) is employed, phase picks will likely be associated in different time windows for a single station. Afterwards, only one(pair) P- and/or S-pick must be selected for depth estimation.

    Args:
         outputs (dict): Dictionary of predicted picks for each station
         outputs_assoc (dict): Dictionary of predicted picks which were associated to the event of interest
         phase_assoc (char): The utilized phase associator. Current options are GaMMA and PyOcto. The latter exports a residual for each associated phase
         mode (char): How a P- or S-phase will be selected. Options are 'max_prob', in which the selected phase will be that at the time window which had the highest probability; and 'min_res', in which the selected phase will be that at the time window which had the lowest residual. This option is only available if the PyOcto phase associator was employed

    Returns:
         outputs_final (Pandas dataframe): A Pandas dataframe with the selected picks for depth estimation 
    """
    statio = []; phas = []; peak_tim = []; peak_valu = []; residua = []; inde = []

    if mode == 'min_res':
       for station in outputs:
           res_p = 9999.; res_s = 9999.
           for secs_bef in outputs_assoc:
               outputs_look = outputs_assoc[secs_bef][0]
               if len(outputs_look[outputs_look['station'] == station]) > 0:
                  outputs_look_phase = outputs_look[outputs_look['station'] == station]
                  if len(outputs_look_phase[outputs_look_phase['phase'] == 'P']) > 0:
                     abs_res_p = abs(outputs_look_phase[outputs_look_phase['phase'] == 'P']['residual'].to_list()[0])
                     if abs_res_p <= abs(res_p):
                        res_p = outputs_look_phase[outputs_look_phase['phase'] == 'P']['residual'].to_list()[0]
                        prob_p = outputs_look_phase[outputs_look_phase['phase'] == 'P']['peak_value'].to_list()[0]
                        t_p = outputs_look_phase[outputs_look_phase['phase'] == 'P']['time'].to_list()[0]
                        pt_p = outputs_look_phase[outputs_look_phase['phase'] == 'P']['peak_time'].to_list()[0]
                        ph_p = outputs_look_phase[outputs_look_phase['phase'] == 'P']['phase'].to_list()[0]
                        in_p = secs_bef
                  if len(outputs_look_phase[outputs_look_phase['phase'] == 'S']) > 0:
                     abs_res_s = abs(outputs_look_phase[outputs_look_phase['phase'] == 'S']['residual'].to_list()[0])
                     if abs_res_s <= abs(res_s):
                        res_s = outputs_look_phase[outputs_look_phase['phase'] == 'S']['residual'].to_list()[0]
                        prob_s = outputs_look_phase[outputs_look_phase['phase'] == 'S']['peak_value'].to_list()[0]
                        t_s = outputs_look_phase[outputs_look_phase['phase'] == 'S']['time'].to_list()[0]
                        pt_s = outputs_look_phase[outputs_look_phase['phase'] == 'S']['peak_time'].to_list()[0]
                        ph_s = outputs_look_phase[outputs_look_phase['phase'] == 'S']['phase'].to_list()[0]
                        in_s = secs_bef
           if abs(res_p) < 9999.:
              statio.append(station)
              phas.append(ph_p)
              peak_tim.append(pt_p)
              peak_valu.append(prob_p)
              residua.append(res_p)
              inde.append(in_p)
           if abs(res_s) < 9999.:
              statio.append(station)
              phas.append(ph_s)
              peak_tim.append(pt_s)
              peak_valu.append(prob_s)
              residua.append(res_s)
              inde.append(in_s)

    if mode == 'max_prob':
       for station in outputs:
           prob_p = 0.0; prob_s = 0.0
           for secs_bef in outputs_assoc:
               outputs_look = outputs_assoc[secs_bef][0]
               if len(outputs_look[outputs_look['station'] == station]) > 0:
                  outputs_look_phase = outputs_look[outputs_look['station'] == station]
                  if len(outputs_look_phase[outputs_look_phase['phase'] == 'P']) > 0:
                     outputs_look_phase_prob_p = outputs_look_phase[outputs_look_phase['phase'] == 'P']
                     if outputs_look_phase_prob_p['peak_value'].to_list()[0] > prob_p:
                        prob_p = outputs_look_phase_prob_p['peak_value'].to_list()[0]
                        t_p = outputs_look_phase_prob_p['time'].to_list()[0]
                        pt_p = outputs_look_phase_prob_p['peak_time'].to_list()[0]
                        ph_p = outputs_look_phase_prob_p['phase'].to_list()[0]
                        if phase_assoc.lower() == 'pyocto':
                           re_p = outputs_look_phase_prob_p['residual'].to_list()[0]
                        in_p = secs_bef
                  if len(outputs_look_phase[outputs_look_phase['phase'] == 'S']) > 0:
                     outputs_look_phase_prob_s = outputs_look_phase[outputs_look_phase['phase'] == 'S']
                     if outputs_look_phase_prob_s['peak_value'].to_list()[0] > prob_s:
                        prob_s = outputs_look_phase_prob_s['peak_value'].to_list()[0]
                        t_s = outputs_look_phase_prob_s['time'].to_list()[0]
                        pt_s = outputs_look_phase_prob_s['peak_time'].to_list()[0]
                        ph_s = outputs_look_phase_prob_s['phase'].to_list()[0]
                        if phase_assoc.lower() == 'pyocto':
                           re_s = outputs_look_phase_prob_s['residual'].to_list()[0]
                        in_s = secs_bef
           if prob_p > 0.0:
              statio.append(station)
              phas.append(ph_p)
              peak_tim.append(pt_p)
              peak_valu.append(prob_p)
              if phase_assoc.lower() == 'pyocto':
                 residua.append(re_p)
              inde.append(in_p)
           if prob_s > 0.0:
              statio.append(station)
              phas.append(ph_s)
              peak_tim.append(pt_s)
              peak_valu.append(prob_s)
              if phase_assoc.lower() == 'pyocto':
                 residua.append(re_s)
              inde.append(in_s)

    if phase_assoc.lower() == 'pyocto':
       df = {'station': statio, 'peak_time': peak_tim, 'peak_value': peak_valu, 'phase': phas, 'residual': residua, 'index': inde}
    else:
         df = {'station': statio, 'peak_time': peak_tim, 'peak_value': peak_valu, 'phase': phas, 'index': inde}
    outputs_final = pd.DataFrame(data=df)

    to_delete = []
    for station in outputs_final['station'].drop_duplicates():
        if len(outputs_final[outputs_final['station'] == station]) == 2:
           out_final_station = outputs_final[outputs_final['station'] == station]
           ts = UTCDateTime(str(out_final_station[out_final_station['phase'] == 'S']['peak_time'].to_list()[0]))
           tp = UTCDateTime(str(out_final_station[out_final_station['phase'] == 'P']['peak_time'].to_list()[0]))
           if ts - tp < 0:
              to_delete.append(station)

    if len(to_delete) > 0:
       for dele in to_delete:
           index = outputs_final[outputs_final['station'] == dele].index
           outputs_final = outputs_final.drop(index=index)

    return outputs_final
