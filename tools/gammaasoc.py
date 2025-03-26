import pandas as pd

from obspy import taup

from gamma.utils import association

from .nicetools import tt_theo_before_assoc
from .utm import from_latlon
from .visualization import plot_assoc


def phase_association(outputs, data, velmod, ev_lon, ev_lat, ev_time, max_dist, plot, mult_windows, secs_before):
    """

    Phase association using GaMMA (https://github.com/AI4EPS/GaMMA). The function receives as input the predicted phase picks obtained from SeisBench (PhaseNet/EQTransformer) and uses a machine-learning model to predict which of them correspond to a certain event.

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
    ev_east, ev_nort, ev_zone_num, ev_zone_lett = from_latlon(ev_lat, ev_lon)

    taupmodel = taup.TauPyModel(model='iasp91')

    # XXX NOTE: Velocity values for config dict below were taken from the BGR crustal velocity model for Germany at 0 km depth, which might not be optimal if we want to automatically associate picks coming from events of diverse depth (earthquakes, explosions, induced). Nevertheless, constant velocities seem to suffice for shallower events.
    config = {'use_amplitude': False, 'dims': ['x(km)', 'y(km)', 'z(km)'], 'vel': {'p': 5.6, 's': 3.9}, 'use_dbscan': False, 'dbscan_min_samples': 10, 'dbscan_eps': 3, 'oversample_factor': 10, 'x(km)': (1e-3 * ev_east - 222.22, 1e-3 * ev_east + 222.22), 'y(km)': (1e-3 * ev_nort - 222.22, 1e-3 * ev_nort + 222.22), 'z(km)': (0, 50), 'max_sigma11': 4.0, 'max_sigma22': 3.0, 'max_sigma12': 3.0, 'method': 'BGMM'}

    config['bfgs_bounds'] = ((config['x(km)'][0] - 1, config['x(km)'][1] + 1), (config['y(km)'][0] - 1, config['y(km)'][1] + 1), (0, config['z(km)'][1] + 1), (None, None))

    if max_dist <= 80.:
       config['min_picks_per_eq'] = 3
       config['min_p_picks_per_eq'] = 2
       config['min_s_picks_per_eq'] = 1
    else:
         config['min_picks_per_eq'] = 6
         config['min_p_picks_per_eq'] = 3
         config['min_s_picks_per_eq'] = 3

    id_picks = []; id_stations = []; station = []; stats = []; phase = []; time = []; peak_value = []; longitude = []; latitude = []; elevation = []; x = []; y = []; z = []
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
               id_picks.append(data[sta]['network'] + '.' + sta)
               station.append(sta)
               phase.append(pick.phase.lower())
               time.append(str(pick.peak_time).replace('T', ' ').replace('Z', ''))
               peak_value.append((pick.peak_value))

               if sta not in stats:
                  id_stations.append(data[sta]['network'] + '.' + sta)
                  stats.append(sta)
                  longitude.append(data[sta]['coords'][1])
                  latitude.append(data[sta]['coords'][0])
                  elevation.append(data[sta]['coords'][2])

                  st_east, st_nort, st_zone_num, st_zone_lett = from_latlon(data[sta]['coords'][0], data[sta]['coords'][1])
                  x.append(st_east * 1e-3)
                  y.append(st_nort * 1e-3)
                  z.append(data[sta]['coords'][2] * 1e-3)

    df = {'id': id_picks, 'station': station, 'timestamp': time, 'prob': peak_value, 'type': phase}
    picks = pd.DataFrame(data=df)

    dfs = {'id': id_stations, 'id_stations': stats, 'longitude': longitude, 'latitude': latitude, 'elevation(m)': elevation, 'x(km)': x, 'y(km)': y, 'z(km)': z}
    stations = pd.DataFrame(data=dfs)

    events, assignments = association(picks, stations, config, method=config['method'])

    events = pd.DataFrame(data=events)
    assignments = pd.DataFrame(data=assignments, columns=['pick_idx', 'event_idx', 'prob_gamma'])

    picks = picks[['timestamp', 'station', 'prob', 'type']].rename(columns={'timestamp': 'peak_time', 'prob': 'peak_value', 'type': 'phase'})
    stations = stations[['id_stations']].rename(columns={'id_stations': 'id'})
    events = events[['event_index', 'time', 'x(km)', 'y(km)', 'z(km)', 'num_picks']].rename(columns={'event_index': 'idx', 'x(km)': 'x', 'y(km)': 'y', 'z(km)': 'z', 'num_picks': 'picks'})

    event_idx_for_assig = []; pick_idx_for_assig = []; station_for_assig = []; phase_for_assig = []; peak_value_for_assig = []; peak_time_for_assig = []
    for evi in range(len(events)):
        for assig in assignments[assignments['event_idx']==evi+1]['pick_idx']:
            event_idx_for_assig.append(evi + 1)
            pick_idx_for_assig.append(assig)
            station_for_assig.append(picks['station'][assig])
            phase_for_assig.append(picks['phase'][assig].upper())
            peak_value_for_assig.append(picks['peak_value'][assig])
            peak_time_for_assig.append(picks['peak_time'][assig])

    assignments = pd.DataFrame(data={'event_idx': event_idx_for_assig, 'pick_idx': pick_idx_for_assig, 'station': station_for_assig, 'phase': phase_for_assig, 'peak_value': peak_value_for_assig, 'peak_time': peak_time_for_assig})

    merged = pd.merge(events, assignments, left_on='idx', right_on='event_idx', suffixes=('', '_pick'))
    merged['idx'] = merged['idx'] - 1

    if plot:
       plot_assoc(ev_time, data, stations, picks, events, merged, mult_windows, secs_before)

    picks_all = len(picks)
    picks_all_p = len(picks[picks['phase'] == 'p'])
    picks_all_s = len(picks[picks['phase'] == 's'])

    print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
    print('GaMMA: For a total of', picks_all, 'predicted picks')
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
        events_assoc[event] = merg[['station', 'time', 'peak_time', 'peak_value', 'phase']]

    if mult_windows:
       return events_assoc, str(secs_before)
    else:
         return events_assoc
