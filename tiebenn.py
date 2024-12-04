import argparse, glob, json, shutil, os
from tools.nonlinloc import inp_files_nlloc_sb, pynlloc, create3dgrid
from tools.retrieve_data import make_station_list
from tools.nicetools import str2bool, strmonth2num
from tools.visualization import plot_picks4loc, plot_hypoc_confidence_ellipsoid
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth

def main(args):

    try:
        event_file = args.event_file
        max_dist = float(args.max_epic_dist)
        picker = args.picker
        nll3d = args.nll3d
        client = args.client
        station_list_db = args.station_database
        min_detections = int(args.min_detections)
        plotpicks = args.plotpicks
        velmod = args.velmod
        denoise = args.denoise
        ph_assoc = args.ph_assoc
        mult_windows = args.mult_windows
        secs_before = int(args.secs_before)

    except:
           class TiebennInputError(Exception):
                 pass
           raise TiebennInputError('Please check the input arguments')

    if picker.lower() not in ['sb_eqt', 'seisbench_eqt', 'seisbench_eqtransformer', 'sb_eqtransformer', 'sb_pn', 'seisbench_pn', 'sb_phasenet', 'seisbench_phasenet']:
       class TiebennPickerError(Exception):
             pass
       raise TiebennPickerError('Incorrect phase-picking model. Accepted inputs are: sb_eqt, sb_eqtransformer, seisbench_eqt, seisbench_eqtransformer, sb_pn, sb_phasenet, seisbench_pn, seisbench_phasenet (not case-sensitive)')
       return

    if client.lower() not in ['sds', 'fdsn']:
       class TiebennClientError(Exception):
             pass
       raise TiebennClientError('Client must be either SDS or FDSN!')
       return

    if min_detections < 3:
       class TiebennMinimumDetectionsError(Exception):
             pass
       raise TiebennMinimumDetectionsError('The event should be detected on at least 3 stations. Please set min_detections accordingly')
       return

    if velmod > 17: # XXX NOTE: This value should change as new velocity models are implemented in the tiebenn/utils directory
       class TiebennVelocityModelError(Exception):
             pass
       raise TiebennVelocityModelError('Selected velocity model is not implemented on Tiebenn')
       return

    if ph_assoc.lower() not in ['gamma', 'pyocto']:
       class PhaseAssociationError(Exception):
             pass
       raise PhaseAssociationError('Unknown phase associator. Options are Gamma and PyOcto')
       return

    if nll3d == True:
       class TiebennNLL3DUnavailable(Exception):
             pass
       raise TiebennNLL3DUnavailable('Depth estimation option using 3D grids is momentarily disabled! Please turn nll3d parameter to False')
       return

    f = open(event_file, 'r')
    for x in f:
        y = x
        x = x.replace('-',' ').replace('_', ' ')
        try:
            ev_year = x.split()[2]
            ev_month = strmonth2num(x.split()[1])
            ev_day = x.split()[0]
            ev_time = ev_year + '-' + ev_month + '-' + ev_day + ' ' + x.split()[3]
        except:
               try:
                   ev_time = y.split()[0] + ' ' + y.split()[1]
               except:
                      class TiebennDatetimeFormatError(Exception):
                            pass
                      raise TiebennDatetimeFormatError('Accepted Datetime formats for the events are: "dd-Mon-yyyy hh:mm:ss" or "yyyy-mm-dd hh:mm:ss"')

        ev_lon = float(x.split()[5])
        ev_lat = float(x.split()[4])

        for rem in glob.glob('*_tiebenn_loc'):
            shutil.rmtree(rem)

        if len(glob.glob('saved_locations')) == 0:
           os.mkdir('saved_locations')

        print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
        print('-----------------------------------------------------------------------')
        print('Detecting and picking P- and S-wave arrival times...')
        print('-----------------------------------------------------------------------')
        print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

        starttime = UTCDateTime(ev_time)

        if mult_windows == False:
           start_t = UTCDateTime(ev_time) - secs_before
           end_t = start_t + 60.
        else:
             if denoise == False:
                secs_before = [10, 8, 6, 4, 2, 0]
             else:
                  secs_before = [10, 5, 2]

        client_list=['BGR', 'LMU', 'GFZ', 'ODC', 'RASPISHAKE', 'RESIF', 'ETH', 'INGV', 'IPGP', 'NIEP', 'ORFEUS']

        try:
           shutil.rmtree(glob.glob('saved_locations/' + str(starttime) + '_tiebenn_loc/')[0])
           print('Location existent in saved_locations directory: the current location will replace the existent one.')
        except:
               pass

        if len(glob.glob(str(starttime) + '_tiebenn_loc/')) > 0:
           shutil.rmtree(glob.glob(str(starttime) + '_tiebenn_loc/')[0])

        print('Making station list...')

        if station_list_db:
           with open('utils/database_clients_stations.json') as data_file:
                data_all = json.load(data_file)

           data = {}
           for cli in client_list:
               for sta in list(data_all[cli].keys()):
                   epicentral_dist = gps2dist_azimuth(data_all[cli][sta]['coords'][0], data_all[cli][sta]['coords'][1], ev_lat, ev_lon)[0] * 0.001
                   if epicentral_dist <= 200.:
                      data[sta] = data_all[cli][sta]

        else:
             stations_json = './tmp_list.json'
             make_station_list(stations_json=stations_json, client_list=client_list, ev_lat=ev_lat, ev_lon=ev_lon, start_time=starttime, end_time=starttime + 60., channel_list=['HH[ZNE12], EH[ZNE12]', 'DN[ZNE12]', 'SH[ZNE12]', 'HG[ZNE12]', 'BH[ZNE12]'], filter_network=[], filter_station=[])

             with open(stations_json) as data_file:
                  data = json.load(data_file)

             os.remove(stations_json)

        if len(data) < min_detections:
           print('STOP: Not enough stations available for P- and S-wave picking within ca 200 km...')
           continue

        print('Removing undesired channels...')
        to_delete = []
        for key in data.keys():
            del_chan = []
            for dele in data[key]['channels']:
                if dele not in ['HHZ', 'HHE', 'HHN', 'EHZ', 'HGZ', 'HGN', 'HGE', 'EHN', 'EHE', 'SHZ', 'SHN', 'SHE', 'DNZ', 'DNE', 'DNN', 'HH1', 'HH2', 'EH1', 'EH2', 'SH1', 'SH2', 'BHZ', 'BHN', 'BHE', 'BH1', 'BH2']:
                   del_chan.append(dele)
            if len(del_chan) > 0:
               for delet in del_chan:
                   data[key]['channels'].remove(delet)
            if len(data[key]['channels']) == 0:
                to_delete.append(key)

        for ht in to_delete:
            data.pop(ht, None)

        if len(data) < min_detections:
           print('STOP: Not enough stations with HH*, EH*, SH*, BH* or DN* channels available within ca 200 km...')
           continue

        print('Calculating epicentral distances...')
        distances_dict = {}
        for et in data:
            distances_dict[et] = "{:.2f}".format(gps2dist_azimuth(data[et]['coords'][0], data[et]['coords'][1], ev_lat, ev_lon)[0] * 0.001)
            data[et]['epic_distance'] = "{:.2f}".format(gps2dist_azimuth(data[et]['coords'][0], data[et]['coords'][1], ev_lat, ev_lon)[0] * 0.001)

        if picker.lower() in ['sb_eqt', 'sb_eqtransformer', 'seisbench_eqt', 'seisbench_eqtransformer', 'sb_pn', 'sb_phasenet', 'seisbench_pn', 'seisbench_phasenet']:
             from tools.sb_tools import picks_sb

             streams = picks_sb(ev_time=ev_time, ev_lon=ev_lon, ev_lat=ev_lat, data=data, max_dist=max_dist, client=client, picker=picker, velmod=velmod, plotpicks=plotpicks, phase_assoc=ph_assoc, pick_sel='min_res', secs_before=secs_before, mult_windows=mult_windows, min_detections=min_detections, denoise=denoise)

        if not glob.glob('*_tiebenn_loc/csv_picks/*.csv'):
           print('Skipping to next event in readfile...')
           continue

        else:
             inp_files_nlloc_sb(ev_lon=ev_lon, ev_lat=ev_lat, ev_time=ev_time, data=data, nll3d=nll3d, velmod=velmod, min_detections=min_detections)

        if not glob.glob(str(starttime) + '_tiebenn_loc/' + 'nlloc_control.in'):
           print('No NLL-control file was produced. Skipping event.')
           try:
               shutil.rmtree(glob.glob('*_tiebenn_loc/')[0])
               pass
           except:
                  pass
           continue

        else:
             control_file = glob.glob('*_tiebenn_loc/')[0] + 'nlloc_control.in'
             control_file_s = glob.glob('*_tiebenn_loc/')[0] + 'nlloc_control_s.in'

             if nll3d:
                if velmod == 6 or velmod == 7:
                   create3dgrid(ev_lon=ev_lon, ev_lat=ev_lat, velmod=velmod)

             sta_gap, sta_nearest = pynlloc(control_file, control_file_s, velmod=velmod, data=data, nll3d=nll3d, plots=plotpicks)

#             if not nll3d:
#                if sta_nearest > 150. or sta_gap > 2000.:
#                   create3dgrid(ev_lon=ev_lon, ev_lat=ev_lat, velmod=velmod)
#                   inp_files_nlloc_sb(ev_lon=ev_lon, ev_lat=ev_lat, ev_time=ev_time, data=data, nll3d=True, velmod=velmod, min_detections=min_detections)
#                   control_file = glob.glob('*_tiebenn_loc/')[0] + 'nlloc_control.in'
#                   control_file_s = glob.glob('*_tiebenn_loc/')[0] + 'nlloc_control_s.in'
#                   sta_gap, sta_nearest = pynlloc(control_file, control_file_s, velmod=velmod, data=data, nll3d=True, plots=plotpicks)

             if plotpicks:
                try:
                    plot_picks4loc(data=data, streams=streams)
                    plot_hypoc_confidence_ellipsoid()
                except:
                       print('Plots unsuccessful. Location rejected?')
                       pass

             shutil.move(glob.glob('*_tiebenn_loc/')[0], 'saved_locations')

    return

def read_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--event_file', type=str, help='The complete path and filename containing the event/s for hypocenter estimation')
    parser.add_argument('--max_epic_dist', type=float, help='The stations used for arrival-time picking will be withing this epicentral distance (in km)')
    parser.add_argument('--picker', default='no', type=str, help='Picker used to phase-picking: SeisBench (eqt or pn)')
    parser.add_argument('--nll3d', default=False, type=str2bool, help='Use 3D velocity model for depth estimation in NonLinLoc')
    parser.add_argument('--client', default='no', type=str, help='FDSN or SDS')
    parser.add_argument('--station_database', type=str2bool, help='True to use json file with stations. False to use ObsPy function get_stations')
    parser.add_argument('--min_detections', default=3, type=int, help='Minimal detections required to hypocenter estimation')
    parser.add_argument('--plotpicks', default=False, type=str2bool, help='If True, it will plot the identified picks on their respective station')
    parser.add_argument('--velmod', type=int, help='Velocity model for phase association and hypocenter estimation with NonLinLoc')
    parser.add_argument('--denoise', default=False, type=str2bool, help='True to use DeepDenoiser on collected waveforms up to a distance of 100 km')
    parser.add_argument('--ph_assoc', default='no', type=str, help='Phase associator. Options are GaMMa and PyOcto (not case sensitive)')
    parser.add_argument('--mult_windows', default=False, type=str2bool, help='Picks in windows with variable start time')
    parser.add_argument('--secs_before', default=0, type=int, help='Seconds before origin time in case mult_window=False')
    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = read_args()
    main(args)
