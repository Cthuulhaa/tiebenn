import glob
import os
import shutil

from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth

from tiebenn.config.input_params import InputParams
from tiebenn.constants.stations_configuration import CHANNELS_TO_KEEP, CLIENT_LIST, CHANNEL_LIST
from tiebenn.logger.logger_settings import logger
from tiebenn.tools.nicetools import calculate_lqs
from tiebenn.tools.nonlinloc import create3dgrid, pynlloc, inp_files_nlloc_sb
from tiebenn.tools.retrieve_data import make_station_list
from tiebenn.tools.sb_tools import picks_sb
from tiebenn.tools.velocity_models import select_velmod
from tiebenn.tools.visualization import plot_hypoc_confidence_ellipsoid, plot_picks4loc, radarplot


def process_event(params: InputParams, ev_lat: float, ev_lon: float, ev_time: str) -> dict:
    if params.vel_mode in ['automatic', 'auto', 'a']:
        params.velmod = select_velmod(ev_lat=ev_lat, ev_lon=ev_lon)

    for rem in glob.glob('*_tiebenn_loc'):
        shutil.rmtree(rem)

    if len(glob.glob('saved_locations')) == 0:
        os.mkdir('saved_locations')

    logger.info('Detecting and picking P- and S-wave arrival times...')

    starttime = UTCDateTime(ev_time)

    if params.mult_windows == False:
        start_t = starttime - params.secs_before
        end_t = start_t + 60.
    else:
        if params.denoise == False:
            params.secs_before = [10, 8, 6, 4, 2, 0]
        else:
            params.secs_before = [10, 5, 2, 0]

    try:
        shutil.rmtree(glob.glob(f"saved_locations/{str(starttime)}_tiebenn_loc/")[0])
        logger.info(
            'Location existent in saved_locations directory: the current location will replace the existent one.')
    except:
        pass

    if len(glob.glob(f"{str(starttime)}_tiebenn_loc/")) > 0:
        shutil.rmtree(glob.glob(f"{str(starttime)}_tiebenn_loc/")[0])

    logger.info('Making station list...')

    data = make_station_list(client_list=CLIENT_LIST, ev_lat=ev_lat, ev_lon=ev_lon, start_time=starttime,
                             end_time=starttime + 60.,
                             channel_list=CHANNEL_LIST, filter_network=[], filter_station=[])

    if len(data) < params.min_detections:
        logger.error('STOP: Not enough stations available for P- and S-wave picking within ca 200 km...')
        return

    logger.info('Removing undesired channels...')

    to_delete = []
    for key in data.keys():
        del_chan = []
        for dele in data[key]['channels']:
            if dele not in CHANNELS_TO_KEEP:
                del_chan.append(dele)
        if len(del_chan) > 0:
            for delet in del_chan:
                data[key]['channels'].remove(delet)
        if len(data[key]['channels']) == 0:
            to_delete.append(key)

    for ht in to_delete:
        data.pop(ht, None)

    if len(data) < params.min_detections:
        logger.error('STOP: Not enough stations with HH*, EH*, SH*, BH* or DN* channels available within ca 200 km...')
        return

    logger.info('Calculating epicentral distances...')
    distances_dict = {}
    for et in data:
        distances_dict[et] = "{:.2f}".format(
            gps2dist_azimuth(data[et]['coords'][0], data[et]['coords'][1], ev_lat, ev_lon)[0] * 0.001)
        data[et]['epic_distance'] = "{:.2f}".format(
            gps2dist_azimuth(data[et]['coords'][0], data[et]['coords'][1], ev_lat, ev_lon)[0] * 0.001)

    if not params.sds_dir:
        logger.info('Reading data using FDSN client.')
        streams = picks_sb(ev_time=ev_time, ev_lon=ev_lon, ev_lat=ev_lat, data=data,
                           max_dist=params.max_epic_dist,
                           client=params.client, picker=params.picker, velmod=params.velmod, plotpicks=params.plots,
                           phase_assoc=params.ph_assoc,
                           pick_sel='max_prob', secs_before=params.secs_before,
                           mult_windows=params.mult_windows,
                           min_detections=params.min_detections, denoise=params.denoise)
    else:
        logger.info("SDS data location: {}. Reading data from sds.".format(params.sds_dir))
        streams = picks_sb(ev_time=ev_time, ev_lon=ev_lon, ev_lat=ev_lat, data=data,
                           max_dist=params.max_epic_dist,
                           client=params.client, picker=params.picker, velmod=params.velmod, plotpicks=params.plots,
                           phase_assoc=params.ph_assoc,
                           pick_sel='max_prob', secs_before=params.secs_before,
                           mult_windows=params.mult_windows,
                           min_detections=params.min_detections, denoise=params.denoise, sds_dir=params.sds_dir)

    if not glob.glob('*_tiebenn_loc/csv_picks/*.csv'):
        logger.info('Skipping to next event in readfile...')
        return

    else:
        inp_files_nlloc_sb(ev_lon=ev_lon, ev_lat=ev_lat, ev_time=ev_time, data=data, nll3d=params.nll3d,
                           velmod=params.velmod,
                           min_detections=params.min_detections)

    if not glob.glob(f"{str(starttime)}_tiebenn_loc/nlloc_control.in"):
        logger.warning('No NLL-control file was produced. Skipping event.')
        try:
            shutil.rmtree(glob.glob('*_tiebenn_loc/')[0])
            pass
        except:
            pass
        return

    else:
        control_file = f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control.in"
        control_file_s = f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control_s.in"

        if params.nll3d:
            if params.velmod == 6 or params.velmod == 7:
                create3dgrid(ev_lon=ev_lon, ev_lat=ev_lat, velmod=params.velmod)

        sta_gap, sta_nearest, event_final_results = pynlloc(control_file, control_file_s, velmod=params.velmod,
                                                            data=data, nll3d=params.nll3d, plots=params.plots)

        if sta_gap == None:
            return

        if params.plots:
            try:
                plot_picks4loc(data=data, streams=streams)
                plot_hypoc_confidence_ellipsoid()
            except:
                logger.warning('Plots unsuccessful. Location rejected?')
                pass
        else:
            for dele in glob.glob('gmt*.gmt'):
                os.remove(dele)

            loc_file = glob.glob('*_tiebenn_loc/loc_eqdatetime*.hyp')[0]
            new_name = f"{glob.glob('*_tiebenn_loc')[0]}/event_location.NLL"
            os.rename(loc_file, new_name)

            for dele in glob.glob('*_tiebenn_loc/loc_eqdatetime*'):
                os.remove(dele)

        lqs_parameters = calculate_lqs(loc_file=f"{glob.glob('*_tiebenn_loc')[0]}/event_location.NLL",
                                       sta_file=f"{glob.glob('*_tiebenn_loc')[0]}/station_coordinates.txt")

        radarplot(lqs_parameters)

        shutil.move(glob.glob('*_tiebenn_loc/')[0], 'saved_locations')

        return event_final_results
