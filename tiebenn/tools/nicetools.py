from collections import defaultdict
import math
import os

import numpy as np
import pandas as pd

from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth


def decimate(phasepicks, data, ev_lon, ev_lat, az_bin, dist_bin, max_bin):
    """

    Detect regions covered by clusters of stations within specific bins and set a maximum amount of stations within each bin, removing the excess.

    Args:
         phasepicks (list): list of phase picks
         data (dict): dictionary that contains---among others---the station coordinates for each station
         ev_lon (float): event's longitude
         ev_lat (float): event's latitude
         az_bin (float): azimuthal bins
         dist_bin (float): epicentral distance bins
         max_bin (int): maximum amount of stations allowed per bin
    Returns:
         filtered_phasepicks (list): list of stations with picks, where there's a maximum of stations per bin
    """
    original_phasepicks = phasepicks.copy()

    if len(original_phasepicks[0]) > 3:
       new_phasepicks, stations_considered = [], []
       for pick in phasepicks:
           if pick[0] not in stations_considered:
              new_phasepicks.append((pick[0], pick[2], pick[3]))
              stations_considered.append(pick[0])

       phasepicks = new_phasepicks

    stations_with_picks = []
    for pick in phasepicks:
        stations_with_picks.append(pick[0])

    distances, azimuths = [], []
    for station in stations_with_picks:
        sta_lat, sta_lon = data[station]['coords'][0], data[station]['coords'][1]
        epic_dist, azimuth, _ = gps2dist_azimuth(ev_lat, ev_lon, sta_lat, sta_lon)

        distances.append(epic_dist * 0.001)
        azimuths.append(azimuth)

    distances = np.array(distances)
    azimuths = np.array(azimuths)

    az_bins = (azimuths // az_bin).astype(int)
    dist_bins = (distances // dist_bin).astype(int)

    bin_keys = list(zip(az_bins, dist_bins))

    binned_stations = defaultdict(list)
    for idx, key in enumerate(bin_keys):
        binned_stations[key].append(idx)

    stations_to_remove = []

    for station_list in binned_stations.values():
        if len(station_list) > max_bin:
           sorted_by_dist = sorted(station_list, key=lambda i: distances[i])
           to_remove = sorted_by_dist[max_bin:]
           stations_to_remove.extend(to_remove)

    stations_to_remove_names = [stations_with_picks[i] for i in stations_to_remove]

    filtered_phasepicks = [pick for pick in original_phasepicks if pick[0] not in stations_to_remove_names]

    return filtered_phasepicks


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


def generate_csv(outputs, data, snr_data, ev_time):
    """

    Generates CSV files with pick detections

    Args:
         outputs (dict): The predicted phases/arrival times obtained from SeisBench for each station
         data (dict): A dictionary with the nearest stations (ObsPy streams) to the event which contain useful channels for phase-picking
         snr_data (dict): Dictionary with the signal-to-noise ratio for each pick at each station
         ev_time (str): ev_time (str): Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
    Returns:
            prediction result files in CSV format for each station
    """
    ev_time = str(UTCDateTime(ev_time))
    os.mkdir(f"{ev_time}_tiebenn_loc/csv_picks")

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
        filename = f"{ev_time}_tiebenn_loc/csv_picks/{stat}_{ev_time}_prediction_results.csv"
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
        string = f"{string}{chanlist[i]},"
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


def station_density(epicenter, stations):
    """

    Calculates station density within a circular area around the epicenter.

    Args:
         epicenter (tuple): (latitude, longitude) for the event
         stations (list): list of tuples [(latitude, longitude), ...] for station locations

    Returns:
         density (float): Station density (stations per square km).
    """
    distances = np.array([gps2dist_azimuth(epicenter[0], epicenter[1], st[0], st[1])[0] / 1000.0
        for st in stations
    ])

    if len(stations) < 5:
       distance = max(distances)
    else:
         distance = np.percentile(distances, 80)

    area = np.pi * (distance ** 2)
    density = len(stations) / area if area > 0 else 0

    return density


def calculate_azimuths(epicenter, stations):
    """

    Calculates azimuths from the epicenter to each station.

    Args:
         epicenter (tuple): (epicenter_latitude, epicenter_longitude)
         stations (list): List of (lat, lon) tuples of the coordinates of each station
    Returns:
         azimuths (NumPy array): Array of station azimuths
    """
    lat_e, lon_e = epicenter

    azimuths = []
    for lat_s, lon_s in stations:
        _, azimuth, _ = gps2dist_azimuth(lat_e, lon_e, lat_s, lon_s)
        azimuths.append(azimuth)

    return np.array(azimuths)


def azimuthal_uniformity_index(azimuths):
    """

    Calculate the Azimuthal Uniformity Index (AUI) for a set of station azimuths.

    Args:
         azimuths (NumPy array): Array of station azimuths

    Returns:
         aui (float): Azimuthal Uniformity Index
    """
    azimuths = np.sort(azimuths)
    gaps = np.diff(np.concatenate(([azimuths[-1] - 360], azimuths)))
    std_gap = np.std(gaps)

    aui = std_gap
    return aui


def azimuthal_gaps(event_lat, event_lon, station_coords):
    """

    Calculates the primary and secondary azimuthal gaps given an event location and station coordinates.

    Args:
         event_lat (float): Epicentral latitude
         event_lon (float): Epicentral longitude
         station_coords (list): List of station coordinates tuples: [(lat1, lon1), (lat2, lon2), (lat3, lon3), ...]

    Returns:
        tuple: (primary_gap, secondary_gap:
    """
    def calculate_azimuth(lat1, lon1, lat2, lon2):
        """
        Calculate the azimuth from one geographic coordinate to another.
        """
        lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
        delta_lon = lon2 - lon1
        x = np.sin(delta_lon) * np.cos(lat2)
        y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon)
        azimuth = np.degrees(np.arctan2(x, y))
        return (azimuth + 360) % 360

    azimuths = [calculate_azimuth(event_lat, event_lon, lat, lon) for lat, lon in station_coords]
    azimuths_sorted = np.sort(azimuths)
    extended_azimuths = np.append(azimuths_sorted, azimuths_sorted[0] + 360)
    gaps = np.diff(extended_azimuths)

    primary_gap = np.max(gaps)
    primary_index = np.argmax(gaps)
    primary_pair = (azimuths_sorted[primary_index], azimuths_sorted[(primary_index + 1) % len(azimuths_sorted)])

    secondary_gap = 0
    secondary_pair = None
    for i in range(len(azimuths_sorted)):
        reduced_azimuths = np.delete(azimuths_sorted, i)
        extended_reduced = np.append(reduced_azimuths, reduced_azimuths[0] + 360)
        reduced_gaps = np.diff(extended_reduced)
        largest_gap = np.max(reduced_gaps)
        if largest_gap > secondary_gap:
            secondary_gap = largest_gap
            secondary_index = np.argmax(reduced_gaps)
            secondary_pair = (reduced_azimuths[secondary_index], reduced_azimuths[(secondary_index + 1) % len(reduced_azimuths)])

    return primary_gap, secondary_gap


def normalize(param, mode, lb, ub):
    """

    Parameter normalization using robust statistics. The resulting normalized parameter ranges between 0 and 1. Each parameter will be normalized after the modes and bundary values of Ramos et al. (2025, in preparation).

    Args:
         param (float): Parameter to be normalized
         mode (str): Choose a normalization mode between 'simple' (simple robust normalization using lb and ub) or 'log' (logarithmic 10-base, robust normalization)
         lb (float): Lower boundary of normalization. If param < lb, then it will be clipped to 0
         ub (float): Upper boundary of normalization. If param > ub, then it will be clipped to 1. Note: lb as well as ub for the parameters to be normalized using the log mode, are actually the logarithm of the normalization boundaries

    Returns:
         norm_param (float): The normalized parameter
    """
    if mode == 'simple':
       norm_param = (np.array(param) - lb) / (ub - lb)
    elif mode == 'log':
       epsilon = 1e-6
       log_param = np.log10(np.array(param) + epsilon)
       norm_param = (log_param - lb) / (ub - lb)

    norm_param = np.clip(norm_param, 0, 1)

    return norm_param


def calculate_lqs(loc_file, sta_file):
    """

    Calculates the Location Quality Score (LQS) metric, as in the work of Ramos et al. (2025, in preparation). This metric describes the quality of automatic locations based on location parameters and station network distribution in the vicinities of the epicenter.

    Args:
         loc_file (str): Full path to location file in NonLinLoc format (see section 'Formats' in http://alomax.free.fr/nlloc/)
         sta_file (str):  Full path to station file in NonLinLoc format (see section 'Formats' in http://alomax.free.fr/nlloc/)

    Returns:
         parameters (Pandas Dataframe): A dataframe with the original event datetime (from the catalog) and the datetime calculated by NonLinLoc, the LQS metric of the event and the eight normalized parameters used for the calculation of the LQS.
    """
    print('Generating radarplot from LQS value...')

    st4loc = []
    with open(loc_file, 'r') as f:
         for line in f:
             if len(line.split()) > 1:
                if line.split()[0] == 'COMMENT':
                   datetime_orig = str(UTCDateTime(f'{line.split()[1]} {line.split()[2]}'.split('"')[1]))
                if line.split()[0] == 'GEOGRAPHIC':
                   li = line.split()
                   epicenter = (float(li[9]), float(li[11]))
                   datetime = UTCDateTime(int(li[2]), int(li[3]), int(li[4]), int(li[5]), int(li[6]), float(li[7]))
                if line.split()[0] == 'QUALITY':
                       rms = float(line.split()[8])
                       npicks = float(line.split()[10])
                       near_sta = float(line.split()[14])
                if line.split()[0] == 'STATISTICS':
                       li = line.split()
                       cov_xx = float(li[8]); cov_xy = float(li[10]); cov_xz = float(li[12])
                       cov_yx = cov_xy; cov_yy = float(li[14]); cov_yz = float(li[16])
                       cov_zx = cov_xz; cov_zy = cov_yz; cov_zz = float(li[18])
                       covariance = np.array([[cov_xx, cov_xy, cov_xz], [cov_yx, cov_yy, cov_yz], [cov_zx, cov_zy, cov_zz]])
                       detcov = np.linalg.det(covariance)
                if line.split()[1] == '?' and float(line.split()[16]) != 0.0:
                   if line.split()[0] not in st4loc:
                      st4loc.append(line.split()[0])
    f.close()

    stations = []
    with open(sta_file, 'r') as f:
         for line in f:
             if line.split()[1] in st4loc:
                st_lat = float(line.split()[3]); st_lon = float(line.split()[4])
                stations.append((st_lat, st_lon))
    f.close()

    density = station_density(epicenter, stations)
    azimuths = calculate_azimuths(epicenter, stations)
    aui = azimuthal_uniformity_index(azimuths)

    azgaps = azimuthal_gaps(epicenter[0], epicenter[1], stations)
    azgap = azgaps[0]
    sec_azgap = azgaps[1]

    norm_detcov = normalize(detcov, 'log', lb=-3.4713284073615576, ub=0.036130328213886044)
    norm_dens = normalize(density, 'log', lb=-3.3078381334234708, ub=-2.4335738944657455)
    norm_azgap = normalize(azgap, 'simple', lb=28.13900043228351, ub=140)
    norm_sazgap = normalize(sec_azgap, 'simple', lb=40.41388529165775, ub=160)
    norm_aui = normalize(aui, 'simple', lb=6.301451360922264, ub=17.577533605788116)
    norm_nsta = normalize(near_sta, 'log', lb=1.000000043429446, ub=1.4771212691961448)

    if rms < 0.128578 and density < 0.0006:
       norm_rms = 1.0
    else:
        norm_rms = normalize(rms, 'log', lb=-0.8781005892308125, ub=-0.2998885611545161)

    norm_npicks = normalize(npicks, 'log', lb=1.6433404945528949, ub=2.0773641073544047)

    weights = {'sta_den': 0.125, 'azgap': 0.05, 'sec_azgap': 0.1, 'rms': 0.1, 'near_sta': 0.1, 'det_cov': 0.35, 'aui': 0.125, 'npicks': 0.05}

    LQS = weights['sta_den'] * norm_dens + weights['azgap'] * (1 - norm_azgap) + weights['sec_azgap'] * (1 - norm_sazgap) + weights['rms'] * (1 - norm_rms) + weights['near_sta'] * (1 - norm_nsta) + weights['det_cov'] * (1 - norm_detcov) + weights['npicks'] * norm_npicks + weights['aui'] * (1 - norm_aui)

    norm_params = {'norm.det.cov': norm_detcov, 'norm.sta.den': norm_dens, 'norm.aui': norm_aui, 'norm.azgap': norm_azgap, 'norm.sec.azgap': norm_sazgap, 'norm.near.sta': norm_nsta, 'norm.rms': norm_rms, 'norm.npicks': norm_npicks}

    merged = {'datetime_orig': datetime_orig, 'datetime': str(datetime), 'LQS': LQS} | norm_params

    parameters = pd.DataFrame([merged])

    return parameters
