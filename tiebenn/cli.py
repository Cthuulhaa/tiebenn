import argparse
import glob
from importlib.resources import files
import os
import shutil

from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth

from tiebenn.tools.nicetools import calculate_lqs, str2bool, strmonth2num
from tiebenn.tools.nonlinloc import create3dgrid, inp_files_nlloc_sb, pynlloc
from tiebenn.tools.retrieve_data import make_station_list
from tiebenn.tools.velocity_models import select_velmod
from tiebenn.tools.visualization import plot_hypoc_confidence_ellipsoid, plot_picks4loc, radarplot


def main(args):

    try:
        event_file = args.event_file
        max_dist = float(args.max_epic_dist)
        picker = args.picker
        nll3d = args.nll3d
        client = args.client
        sds_dir = args.sds_dir
        min_detections = round(args.min_detections)
        plots = args.plots
        vel_mode = args.vel_mode
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

    if client.lower() == 'sds':
       if not sds_dir:
          class TiebennSDSDirError(Exception):
               pass
          raise TiebennSDSDirError('Parameter sds_dir must be set')
          return

    if min_detections < 3:
       print('WARNING! Minimal detections too small. Set to be 3.')
       min_detections = 3

    if vel_mode.lower() in ['manual', 'man', 'm']:
       if velmod is None:
          class TiebennVelocityModeError(Exception):
                pass
          raise TiebennVelocityModeError('Velocity model must be selected when manual velocity mode in use.')
          return

       if velmod not in [6, 7, 12, 13, 17]:
          velmod_name = f"v{str(velmod)}"

          model_path = files('tiebenn.data.velocity_models').joinpath(velmod_name)

          if model_path.is_file():
             print(f"Velocity model {str(velmod)} selected")
          else:
               class VelModLoadError(Exception):
                     pass
               raise VelModLoadError(f"Velocity model {str(velmod)} does not exist.")
       elif velmod in [6, 7]:
            print('Crust1.0 model will be used for seismic location.')
       else:
            class TiebennNLL3DUnavailable(Exception):
                  pass
            raise TiebennNLL3DUnavailable('Depth estimation option using 3D grids is momentarily disabled!')
            return
    elif vel_mode.lower() not in ['manual', 'man', 'm', 'automatic', 'auto', 'a']:
         class TiebennVelocityModeError(Exception):
               pass
         raise TiebennVelocityModeError('vel_mode must be set to either manual (m) or automatic (a).')
         return

    if ph_assoc.lower() not in ['gamma', 'g', 'pyocto', 'p']:
       class PhaseAssociationError(Exception):
             pass
       raise PhaseAssociationError('Unknown phase associator. Options are Gamma (g) and PyOcto (p)')
       return

    if nll3d == True:
       class TiebennNLL3DUnavailable(Exception):
             pass
       raise TiebennNLL3DUnavailable('Depth estimation option using 3D grids is momentarily disabled! Please turn nll3d parameter to False')
       return

    f = open(event_file, 'r')
    for x in f:
        x = x.replace('T', ' ')
        y = x
        x = x.replace('-',' ').replace('_', ' ')
        try:
            ev_year = x.split()[2]
            ev_month = strmonth2num(x.split()[1])
            ev_day = x.split()[0]
            ev_time = f"{ev_year}-{ev_month}-{ev_day} {x.split()[3]}"
        except:
               try:
                   ev_time = f"{y.split()[0]} {y.split()[1]}"
               except:
                      class TiebennDatetimeFormatError(Exception):
                            pass
                      raise TiebennDatetimeFormatError('Accepted Datetime formats for the events are: "dd-Mon-yyyy hh:mm:ss", "yyyy-mm-ddThh:mm:ss", or "yyyy-mm-dd hh:mm:ss"')

        ev_lon = float(x.split()[5])
        ev_lat = float(x.split()[4])

        if vel_mode in ['automatic', 'auto', 'a']:
           velmod = select_velmod(ev_lat=ev_lat, ev_lon=ev_lon)

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
                  secs_before = [10, 5, 2, 0]

        client_list=['BGR', 'LMU', 'GFZ', 'ODC', 'RASPISHAKE', 'RESIF', 'ETH', 'INGV', 'IPGP', 'NIEP', 'ORFEUS']

        try:
           shutil.rmtree(glob.glob(f"saved_locations/{str(starttime)}_tiebenn_loc/")[0])
           print('Location existent in saved_locations directory: the current location will replace the existent one.')
        except:
               pass

        if len(glob.glob(f"{str(starttime)}_tiebenn_loc/")) > 0:
           shutil.rmtree(glob.glob(f"{str(starttime)}_tiebenn_loc/")[0])

        print('Making station list...')

        data = make_station_list(client_list=client_list, ev_lat=ev_lat, ev_lon=ev_lon, start_time=starttime, end_time=starttime + 60., channel_list=['HH[ZNE12], EH[ZNE12]', 'DN[ZNE12]', 'SH[ZNE12]', 'HG[ZNE12]', 'BH[ZNE12]'], filter_network=[], filter_station=[])

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
           from tiebenn.tools.sb_tools import picks_sb

           if not sds_dir:
              streams = picks_sb(ev_time=ev_time, ev_lon=ev_lon, ev_lat=ev_lat, data=data, max_dist=max_dist, client=client, picker=picker, velmod=velmod, plotpicks=plots, phase_assoc=ph_assoc, pick_sel='max_prob', secs_before=secs_before, mult_windows=mult_windows, min_detections=min_detections, denoise=denoise)
           else:
                streams = picks_sb(ev_time=ev_time, ev_lon=ev_lon, ev_lat=ev_lat, data=data, max_dist=max_dist, client=client, picker=picker, velmod=velmod, plotpicks=plots, phase_assoc=ph_assoc, pick_sel='max_prob', secs_before=secs_before, mult_windows=mult_windows, min_detections=min_detections, denoise=denoise, sds_dir=sds_dir)

        if not glob.glob('*_tiebenn_loc/csv_picks/*.csv'):
           print('Skipping to next event in readfile...')
           continue

        else:
             inp_files_nlloc_sb(ev_lon=ev_lon, ev_lat=ev_lat, ev_time=ev_time, data=data, nll3d=nll3d, velmod=velmod, min_detections=min_detections)

        if not glob.glob(f"{str(starttime)}_tiebenn_loc/nlloc_control.in"):
           print('No NLL-control file was produced. Skipping event.')
           try:
               shutil.rmtree(glob.glob('*_tiebenn_loc/')[0])
           except:
                  pass
           continue

        else:
             control_file = f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control.in"
             control_file_s = f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control_s.in"

             if nll3d:
                if velmod == 6 or velmod == 7:
                   create3dgrid(ev_lon=ev_lon, ev_lat=ev_lat, velmod=velmod)

             sta_gap, sta_nearest = pynlloc(control_file, control_file_s, velmod=velmod, data=data, nll3d=nll3d, plots=plots)

             if sta_gap == None:
                continue

#             if not nll3d:
#                if sta_nearest > 150. or sta_gap > 2000.:
#                   create3dgrid(ev_lon=ev_lon, ev_lat=ev_lat, velmod=velmod)
#                   inp_files_nlloc_sb(ev_lon=ev_lon, ev_lat=ev_lat, ev_time=ev_time, data=data, nll3d=True, velmod=velmod, min_detections=min_detections)
#                   control_file = glob.glob('*_tiebenn_loc/')[0] + 'nlloc_control.in'
#                   control_file_s = glob.glob('*_tiebenn_loc/')[0] + 'nlloc_control_s.in'
#                   sta_gap, sta_nearest = pynlloc(control_file, control_file_s, velmod=velmod, data=data, nll3d=True, plots=plots)

             if plots:
                try:
                    plot_picks4loc(data=data, streams=streams)
                    plot_hypoc_confidence_ellipsoid()
                except:
                       print('Plots unsuccessful. Location rejected?')
                       pass
             else:
                  for dele in glob.glob('gmt*.gmt'):
                      os.remove(dele)

                  loc_file = glob.glob('*_tiebenn_loc/loc_eqdatetime*.hyp')[0]
                  new_name = f"{glob.glob('*_tiebenn_loc')[0]}/event_location.NLL"
                  os.rename(loc_file, new_name)

                  for dele in glob.glob('*_tiebenn_loc/loc_eqdatetime*'):
                      os.remove(dele)

             lqs_parameters = calculate_lqs(loc_file=f"{glob.glob('*_tiebenn_loc')[0]}/event_location.NLL", sta_file=f"{glob.glob('*_tiebenn_loc')[0]}/station_coordinates.txt")

             radarplot(lqs_parameters)

             shutil.move(glob.glob('*_tiebenn_loc/')[0], 'saved_locations')


    return

def read_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--event_file', type=str, help='The complete path and filename containing the event/s for hypocenter estimation')
    parser.add_argument('--max_epic_dist', type=float, help='The stations used for arrival-time picking will be withing this epicentral distance (in km)')
    parser.add_argument('--picker', default='PhaseNet', type=str, help='Picker used to phase-picking: SeisBench (eqt or pn)')
    parser.add_argument('--nll3d', default=False, type=str2bool, help='Use 3D velocity model for depth estimation in NonLinLoc')
    parser.add_argument('--client', default='FDSN', type=str, help='FDSN or SDS')
    parser.add_argument('--sds_dir', default='/', type=str, help='full path to SeisComp3 directory')
    parser.add_argument('--min_detections', default=3, type=int, help='Minimal detections required to hypocenter estimation')
    parser.add_argument('--plots', default=False, type=str2bool, help='If True, it will plot several figures')
    parser.add_argument('--vel_mode', default='auto', type=str, help='Velocity model mode: automatic or manual selection. If manual, velmod must be specified')
    parser.add_argument('--velmod', type=int, help='Velocity model for phase association and hypocenter estimation with NonLinLoc')
    parser.add_argument('--denoise', default=False, type=str2bool, help='True to use DeepDenoiser on collected waveforms up to a distance of 100 km')
    parser.add_argument('--ph_assoc', default='no', type=str, help='Phase associator. Options are GaMMa and PyOcto (not case sensitive)')
    parser.add_argument('--mult_windows', default=False, type=str2bool, help='Picks in windows with variable start time')
    parser.add_argument('--secs_before', default=0, type=int, help='Seconds before origin time in case mult_window=False')
    args = parser.parse_args()

    return args


def main_cli():

    args = read_args()
    main(args)
