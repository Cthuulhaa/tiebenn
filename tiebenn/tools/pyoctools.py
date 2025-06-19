import datetime
import os

import numpy as np
import pandas as pd
import pyocto

from obspy import UTCDateTime, taup

from .nicetools import tt_theo_before_assoc
from .velocity_models import velmods
from .visualization import plot_assoc


def phase_association(outputs, data, velmod, ev_lon, ev_lat, ev_time, max_dist, plot, mult_windows, secs_before):
    """

    Phase association using PyOcto (https://github.com/yetinam/pyocto). The function receives as input the predicted phase picks from SeisBench (PhaseNet/EQTransformer), uses 4D space-time partitioning and can employ homogeneous and 1D velocity models.

    NOTE: A future development for this function could be the case where the phase association yields more than one event, not originally detected by the AD-Detector. In this case, additional outputs would be required to account for eventual new events

    Args:
         outputs (dict): The predicted phases/arrival times obtained from SeisBench for each station
         data (dict): A dictionary with information regarding the stations on which picks were predicted
         velmod (int): The velocity model used for phase association
         ev_lon (float): Longitude of the located event
         ev_lat (float): Latitude of the located event
         ev_time (str) Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
         plot (bool): If true, it will plot the detected event(s) and the corresponding associated picks
         mult_windows (bool): If picks were predicted in the multiple-windows-mode, then the code will name de output plots accordingly
         secs_before (int): Seconds before the event time to start retrieving waveforms, among other uses
    Returns:
            events_assoc (dict): A dictionary with all the predicted events and their associated phases in a pandas dataframe
            PhAssoc_event<event_number>.pdf: A visualization of the associated picks

    """
    taupmodel = taup.TauPyModel(model='iasp91')

    station = []; phase = []; time = []; stats = []; longitude = []; latitude = []; elevation = []; peak_value = []

    for sta in outputs:
        arrivals = taupmodel.get_travel_times(source_depth_in_km=5., distance_in_degree=float(data[sta]['epic_distance'])/111.11)
        for i in arrivals:
            if i.name.lower() == 's':
               s_secs_taup = float(i.time)
            if i.name.lower() == 'p':
               p_secs_taup = float(i.time)

        for pick in outputs[sta].picks:
            teo_verdict = tt_theo_before_assoc(ev_time=ev_time, teo_p_time=p_secs_taup, teo_s_time=s_secs_taup, pick=pick, tol_p=6, tol_s=7) #XXX NOTE: The tolerance indicates that if for a pick the theoretical time (for an event at 5 km depth) differs <tolerance> seconds from the observed time, the observation will be discarded. This aims to eliminate false positives, but it is problematic e.g. for small events or low SNR, since more false positives are expected and the association might fail with too few picks. For this reason, the tolerance must allow an error margin in the picks

            if teo_verdict == 'pass':
               station.append(pick.trace_id.split('.')[1])
               phase.append(pick.phase)
               time.append(np.datetime64(pick.peak_time))
               peak_value.append(pick.peak_value)

               if sta not in stats:
                  stats.append(sta)
                  longitude.append(data[sta]['coords'][1])
                  latitude.append(data[sta]['coords'][0])
                  elevation.append(data[sta]['coords'][2])

    if len(station) <= 3:
       print('PyOcto: not enough picks passed the comparison with theoretical arrival times! Skipping.')
       return

    stations = pd.DataFrame(data={'id': stats, 'longitude': longitude, 'latitude': latitude, 'elevation': elevation})
    picks = pd.DataFrame(data={'station': station, 'phase': phase, 'time': time, 'peak_value': peak_value})

    picks['peak_time'] = picks['time']
    picks['time'] = picks['time'].apply(lambda x: x.timestamp())

    if velmod in [12, 13, 17]:
       velmod_ = 7
       layers_from_velmod = velmods(model=velmod_, ev_lon=ev_lon, ev_lat=ev_lat)

    else:
         layers_from_velmod = velmods(model=velmod, ev_lon=ev_lon, ev_lat=ev_lat)

    depth = []; vp = []; vs = []
    for layer in layers_from_velmod:
        depth.append(float(layer.split()[1]))
        vp.append(float(layer.split()[2]))
        vs.append(float(layer.split()[4]))

    layers = pd.DataFrame(data={'depth': depth, 'vp': vp, 'vs': vs}).sort_values('depth')
    velmod_path = f"velocity_model{str(secs_before)}"
    pyocto.VelocityModel1D.create_model(layers, 1., 300, 300, velmod_path)
    velocity_model = pyocto.VelocityModel1D(velmod_path, tolerance=2.5, association_cutoff_distance=300.)

    outputs_p = 0; outputs_s = 0
    for picks_ in outputs:
        for phase in outputs[picks_].picks:

            if phase.phase.lower() == 'p':
               outputs_p = outputs_p + 1
            else:
                 outputs_s = outputs_s + 1

    if outputs_s <= 3:
       associator = pyocto.OctoAssociator.from_area(lat=(ev_lat - 1.3, ev_lat + 1.3), lon=(ev_lon - 1.3, ev_lon + 1.3), zlim=(0, 100), time_before=300, velocity_model=velocity_model, n_p_picks=2, n_s_picks=2, n_picks=3, n_p_and_s_picks=1, pick_match_tolerance=2.5)
    elif outputs_s > 3 and outputs_s < 6:
       associator = pyocto.OctoAssociator.from_area(lat=(ev_lat - 1.3, ev_lat + 1.3), lon=(ev_lon - 1.3, ev_lon + 1.3), zlim=(0, 100), time_before=400, velocity_model=velocity_model, n_p_picks=2, n_s_picks=2, n_picks=5, n_p_and_s_picks=2, pick_match_tolerance=2.5)
    else:
         associator = pyocto.OctoAssociator.from_area(lat=(ev_lat - 1.3, ev_lat + 1.3), lon=(ev_lon - 1.3, ev_lon + 1.3), zlim=(0, 100), time_before=400, velocity_model=velocity_model, n_p_picks=3, n_s_picks=3, n_picks=6, n_p_and_s_picks=3, pick_match_tolerance=2.5)

    associator.transform_stations(stations)
    events, assignments = associator.associate(picks, stations)
    associator.transform_events(events)

    try:
        events['time'] = events['time'].apply(datetime.datetime.fromtimestamp, tz=datetime.timezone.utc)
        merged = pd.merge(events, assignments, left_on='idx', right_on='event_idx', suffixes=('', '_pick'))
    except:
           print('PyOcto: No arrivals were associated to this event %s seconds before the event.' % secs_before)
           return

    if plot:
       plot_assoc(ev_time, data, stations, picks, events, merged, mult_windows, secs_before)

    picks_all = len(picks)
    picks_all_p = len(picks[picks['phase'] == 'P'])
    picks_all_s = len(picks[picks['phase'] == 'S'])

    print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
    print('PyOcto: Phase association with waveforms starting %s seconds before event time' % str(secs_before))
    print('------------------------------------------------------------------------------')
    print('For a total of %s predicted picks which are comparable to theoretical picks' % str(picks_all))
    print('of which', picks_all_p, 'are P-picks and', picks_all_s, 'are S-picks:')
    print('%', len(events), 'event(s) was/were detected')

    for event in range(len(events)):
        print('% for event', event + 1, ':')
        for_event = merged[merged['idx'] == event]
        picks_event = len(for_event)
        picks_event_p = len(for_event[for_event['phase'] == 'P'])
        picks_event_s = len(for_event[for_event['phase'] == 'S'])
        print('%', picks_event, 'picks were associated')
        print('% of which', picks_event_p, 'are P-picks and', picks_event_s, 'are S-picks')
    print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

    events_assoc = {}
    for event in range(len(events)):
        merg = merged[merged['idx'] == event]
        events_assoc[event] = merg[['station', 'time', 'peak_time', 'peak_value', 'phase', 'residual']]

    os.remove(velmod_path)

    if mult_windows:
       return events_assoc, str(secs_before)
    else:
         return events_assoc
