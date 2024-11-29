import glob, os
from obspy import UTCDateTime, read
import pygmt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plotpicks_sb_mw(data, streams, starttime, predictions, picks, picks_final):
    """

    Generates a plot which includes time windows showing the waveforms on the different channels on every station used for wave detection under the multiple time-window approach. The predicted P- and/or S-picks are shown as thicks: the ones used for the depth estimation are colored and the rest -- obtained from predictions with other time offsets (and/or not associated to the event) -- are shown thinner and in grey. The probability function of the P- or S-picks is also plotted in colored lines, while that obtained from predictions with other time offsets are in grey. This function is adapted to the outputs produced by the pickers in SeisBench. 

    Args:
         data (dict): A dictionary with station information
         streams (dict): A dictionary with the ObsPy streams for each station
         starttime (str): Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
         predictions (dict): A dictionary with the probability function for a detection, P- or S-wave on a given station, as obtained from the Annotate method implemented on SeisBench
         picks (dict): A dictionary with information and classification of P- or S-picks on each station, as obtained from the Classify method implemented on SeisBench
         picks_final (dict):

    Returns:
         Figures in PDF format, which are stored in the saved_locations/<ev_time>_tiebenn_loc directory
    """
    secs_bef = 10

    for sta in picks:
        try:
            ax0 = 0
            fig, axis = plt.subplots(len(streams[sta]) + 1, 1, figsize=(15, 7), sharex=True, tight_layout=True)
            for ch in range(len(streams[sta])):
                axis[ax0].plot(streams[sta][ch].times() - secs_bef, streams[sta][ch].data / np.amax(streams[sta][ch].data), 'k', linewidth=1)
                axis[ax0].set_ylabel(streams[sta][ch].stats.channel)

                for seconds in picks[sta]:
                    for pick in picks[sta][seconds].picks:
                        pick_time = pick.peak_time - streams[sta][ch].stats.starttime - secs_bef
                        axis[ax0].vlines(pick_time, min(streams[sta][ch].data / np.amax(streams[sta][ch].data)), max(streams[sta][ch].data / np.amax(streams[sta][ch].data)), linestyles='solid', colors='lightgrey', linewidth=1)

                pick_pl = picks_final[picks_final['station'] == sta]
                for pick in range(len(pick_pl)):
                    pick_time = UTCDateTime(str(pick_pl['peak_time'].to_list()[pick])) - streams[sta][ch].stats.starttime - secs_bef
                    if pick_pl['phase'].to_list()[pick] == 'P':
                       axis[ax0].vlines(pick_time, min(streams[sta][ch].data / np.amax(streams[sta][ch].data)), max(streams[sta][ch].data / np.amax(streams[sta][ch].data)), linestyles ='solid', colors ='c', label='P-pick')
                    if pick_pl['phase'].to_list()[pick] == 'S':
                       axis[ax0].vlines(pick_time, min(streams[sta][ch].data / np.amax(streams[sta][ch].data)), max(streams[sta][ch].data / np.amax(streams[sta][ch].data)), linestyles ='solid', colors ='r', label='S-pick')
                    axis[ax0].legend(loc='upper right')

                if ax0 == 0:
                   axis[ax0].set_title(streams[sta][ch].stats.network + '.' + sta + ' ' + str(starttime))
                ax0 = ax0 + 1

            for seconds in predictions[sta]:
                for preds in range(len(predictions[sta][seconds])):
                    if predictions[sta][seconds][preds].stats.channel.split('_')[-1] != 'N':
                       offset = predictions[sta][seconds][preds].stats.starttime - streams[sta][ch].stats.starttime
                       axis[ax0].plot(predictions[sta][seconds][preds].times() + offset - secs_bef, predictions[sta][seconds][preds].data, color='lightgrey', linewidth=1)

            for pick in range(len(pick_pl)):
                ind_max = pick_pl['index'].to_list()[pick]
                for pred in range(len(predictions[sta][ind_max])):
                    if predictions[sta][str(ind_max)][pred].stats.channel.split('_')[-1] == pick_pl['phase'].to_list()[pick]:
                       offset = predictions[sta][str(ind_max)][pred].stats.starttime - streams[sta][ch].stats.starttime
                       axis[ax0].plot(predictions[sta][str(ind_max)][pred].times() + offset - secs_bef, predictions[sta][str(ind_max)][pred].data, label=predictions[sta][str(ind_max)][pred].stats.channel.split('_')[-1])
            axis[ax0].text(0, -0.3, 'Epicentral distance: ' + data[sta]['epic_distance'] + ' (km)', fontsize=12)
            axis[ax0].legend()
            axis[ax0].set_ylim(0, 1)
            axis[ax0].set_ylabel('Probability')
            axis[ax0].set_xlabel('Seconds from start time')
        except:
               print('No data to plot for station', sta)
               pass
        try:
            os.mkdir(str(starttime) + '_tiebenn_loc/plots_picks')
        except:
               pass
        try:
            savename = str(starttime) + '_tiebenn_loc/plots_picks/' + streams[sta][0].stats.network + '.' + sta + str(starttime) + '.pdf'
            plt.savefig(savename)
            plt.close()
        except:
               pass

    return


def plotpicks_sb(data, streams, starttime, predictions, picks):
    """

    Generates a plot which includes time windows showing the waveforms on the different channels on every station used for wave detection. The identified P- and/or S-picks are shown as thicks. The probability function of the detected event and P- or S-picks is also plotted. This function is adapted to the outputs produced by the pickers in SeisBench. 

    Args:
         data (dict): A dictionary with station information
         streams (dict): A dictionary with the ObsPy streams for each station
         starttime (str): Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
         predictions (dict): A dictionary with the probability function for a detection, P- or S-wave on a given station, as obtained from the Annotate method implemented on SeisBench
         picks (dict): A dictionary with information and classification of P- or S-picks on each station, as obtained from the Classify method implemented on SeisBench

    Returns:
         Figures in PDF format, which are stored in the saved_locations/<ev_time>_tiebenn_loc directory
    """

    for sta in streams:
        try:
            ax0 = 0
            fig, axis = plt.subplots(len(streams[sta]) + 1, 1, figsize=(15, 7), sharex=True, tight_layout=True)
            for ch in range(len(streams[sta])):
                axis[ax0].plot(streams[sta][ch].times(), streams[sta][ch].data / np.amax(streams[sta][ch].data), 'k', linewidth=1)
                for pick in picks[sta].picks:
                    if pick.phase == 'P':
                       axis[ax0].vlines(pick.peak_time - streams[sta][ch].stats.starttime, min(streams[sta][ch].data / np.amax(streams[sta][ch].data)), max(streams[sta][ch].data / np.amax(streams[sta][ch].data)), linestyles='solid', colors='c', label='P-pick')
                       axis[ax0].legend(loc='upper right')
                    if pick.phase == 'S':
                       axis[ax0].vlines(pick.peak_time - streams[sta][ch].stats.starttime, min(streams[sta][ch].data / np.amax(streams[sta][ch].data)), max(streams[sta][ch].data / np.amax(streams[sta][ch].data)), linestyles='solid', colors ='r', label='S-pick')
                       axis[ax0].legend(loc='upper right')
                axis[ax0].set_ylabel(streams[sta][ch].stats.channel)
                if ax0 == 0:
                   axis[ax0].set_title(streams[sta][ch].stats.network + '.' + sta + ' ' + str(starttime))
                ax0 = ax0 + 1
            for preds in range(len(predictions[sta])):
                if predictions[sta][preds].stats.channel.split('_')[-1] != 'N':
                   offset = predictions[sta][preds].stats.starttime - streams[sta][ch].stats.starttime
                   axis[ax0].plot(predictions[sta][preds].times() + offset, predictions[sta][preds].data, label=predictions[sta][preds].stats.channel.split('_')[-1])
            axis[ax0].text(0, -0.3, 'Epicentral distance: ' + data[sta]['epic_distance'] + ' (km)', fontsize=12)
            axis[ax0].legend()
            axis[ax0].set_ylim(0, 1)
            axis[ax0].set_ylabel('Probability')
            axis[ax0].set_xlabel('Seconds from start time')
        except:
               print('No data to plot for station', sta)
               pass
        try:
            os.mkdir(str(starttime) + '_tiebenn_loc/plot_picks')
        except:
               pass
        try:
            savename = str(starttime) + '_tiebenn_loc/plot_picks/' + streams[sta][0].stats.network + '.' + sta + str(starttime) + '.pdf'
            plt.savefig(savename)
            plt.close()
        except:
               pass

    return


def plot_assoc(ev_time, data, stations, picks, events, merged, mult_windows, secs_before):
    """

    The visualization of the associated picks using PyOcto

    Args:
         ev_time (str): Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
         data (dict): A dictionary with information about the stations with predicted picks
         stations (pandas dataframe): A dataframe with the coordinates and elevation for each station
         picks (pandas dataframe): Information for each predicted pick: station, phase (P or S), time and probability
         events (pandas dataframe): The event(s) to which PyOcto associated the picks
         merged (pandas dataframe): A single dataframe containing all the information for the associated picks
         mult_windows (bool): If picks were predicted in the multiple-windows-mode, then the code will name de output plots accordingly
         secs_before (int): Seconds before the event time to start retrieving waveforms, among other uses

    Returns:
            PhAssoc_event<event_number>.pdf: Figure in PDF format, which is stored in the saved_locations/<ev_time>_tiebenn_loc directory
    """
    for ev in range(len(events)):
        max_dist = 0; max_time = 0
        for pi in range(len(picks)):
            t_pick = UTCDateTime(str(picks['peak_time'][pi]))
            time = t_pick - ev_time
            if time >= max_time:
               max_time = time
            st = picks['station'][pi]
            epic_dist = float(data[st]['epic_distance'])
            if epic_dist >= max_dist:
               max_dist = epic_dist
            if picks['phase'][pi].lower() == 'p':
               plt.plot(epic_dist, time, 'c', marker='x', markersize=6)
            if picks['phase'][pi].lower() == 's':
               plt.plot(epic_dist, time, 'r', marker='x', markersize=6)

        merged_plot = merged[merged['idx'] == ev]
        merged_indexes = merged.index[merged['idx'] == ev].to_list()

        for merg in merged_indexes:
            st = merged_plot['station'][merg]
            epic_dist = float(data[st]['epic_distance'])
            time = UTCDateTime(str(merged_plot['peak_time'][merg])) - ev_time
            if merged_plot['phase'][merg].lower() == 'p':
               plt.plot(epic_dist, time, 'c', marker='.', markersize=25, markeredgewidth=1.5, markeredgecolor='k')
            if merged_plot['phase'][merg].lower() == 's':
               plt.plot(epic_dist, time, 'r', marker='.', markersize=25, markeredgewidth=1.5, markeredgecolor='k')

        for station in stations['id']:
            epic_dist = float(data[station]['epic_distance'])
            plt.text(epic_dist, 0, station, size=10, rotation=40, ha='center', va='center', bbox=dict(boxstyle='square', ec=(0, 0, 0), fc=(1, 0.9, 0.9)))

        plt.plot(-10, -10, 'cx', markersize=6, label='Non-associated P')
        plt.plot(-10, -10, 'rx', markersize=6, label='Non-associated S')
        plt.plot(-10, -10, 'c.', markersize=25, markeredgewidth=1.5, markeredgecolor='k', label='Associated P')
        plt.plot(-10, -10, 'r.', markersize=25, markeredgewidth=1.5, markeredgecolor='k', label='Associated S')
        plt.xlabel('Epicentral distance [km]')
        plt.ylabel('Seconds after event')
        plt.xlim(0, max_dist + 5)
        plt.ylim(0, max_time + 5)
        plt.grid(True)
        plt.legend(loc='best')
        title = 'Picks associated to event ' + str(ev + 1)
        plt.title(title)
        if mult_windows == True:
           savename = 'PhAssoc_event' + str(ev + 1) + '_' + str(secs_before) + '.pdf'
        else:
             savename = 'PhAssoc_event' + str(ev + 1) + '.pdf'
        plt.savefig(savename)
        plt.close()


def epic_sta_plot():
    """

    Generates a plot of the event's epicenter and the stations with picks used for depth estimation with NonLinLoc.
    """
    stations_file = glob.glob('*_tiebenn_loc/')[0] + 'station_coordinates.txt'

    locfile = glob.glob('*_tiebenn_loc/loc_eqdatetime.*.*.grid0.loc.hyp')[0]

    sta_list = {}
    with open(stations_file) as stats:
         for line in stats:
             sta_list[line.split()[1]] = [float(line.split()[3]), float(line.split()[4])]

    station = []; latitude = []; longitude = []
    with open(locfile) as locf:
         for line in locf:
             if len(line.split()) != 0:
                if line.split()[0] == 'GEOGRAPHIC':
                   ev_longitude = float(line.split()[11])
                   ev_latitude = float(line.split()[9])
                if line.split()[0] in sta_list and float(line.split()[16]) != 0.0:
                   station.append(line.split()[0])
                   latitude.append(sta_list[line.split()[0]][0])
                   longitude.append(sta_list[line.split()[0]][1])

    df = {'station': station, 'latitude': latitude, 'longitude': longitude}
    stations = pd.DataFrame(data=df)

    fig = pygmt.Figure()
    fig.coast(region=[ev_longitude - 1.3, ev_longitude + 1.3, ev_latitude - 1.3, ev_latitude + 1.3], land='lightgray', water='white', borders='1/1p', frame='ag')
    fig.plot(x=ev_longitude, y=ev_latitude, style='s0.3c', fill='darkred', pen='black', label='Epicenter')
    fig.text(x=stations.longitude, y=stations.latitude, text=stations.station, font='12p,Courier-Bold,black', justify='LT')
    fig.legend()

    fig.savefig(glob.glob('*_tiebenn_loc/')[0] + 'epicenter_stations.pdf')

    return


def plot_picks4loc(data, streams):
    """
    Generates a plot depicting all the stations sorted by epicentral distance with picks used for depth estimation with NonLinLoc
    """
    stations_file = glob.glob('*_tiebenn_loc/')[0] + 'station_coordinates.txt'
    observations_file = glob.glob('*_tiebenn_loc/')[0] + 'obs_ttimes.obs'
    secs_bef = 10

    nll_loctime = glob.glob('*_tiebenn_loc/*.hyp')[0]
    event_time_nll = UTCDateTime(nll_loctime.split('.')[2] + nll_loctime.split('.')[3])

    fig, ax = plt.subplots(1, 1, figsize=(11, 15), tight_layout=True)
    for st in open(stations_file, 'r'):
        station = st.split()[1]
        distance = float(data[station]['epic_distance'])
        waveform = streams[station][0]
        time_diff = streams[station][0].stats.starttime + secs_bef - event_time_nll

        ax.plot(streams[station][0].times() - time_diff, (streams[station][0].data / np.amax(streams[station][0].data)) + distance, 'k', linewidth=1)
        ax.text(0, distance + 0.1, station)

        for obs in open(observations_file, 'r'):
            if obs.split()[0] == station:
               year = int(obs.split()[6][0:4])
               month = int(obs.split()[6][4:6])
               day = int(obs.split()[6][6:8])
               hour = int(obs.split()[7][0:2])
               minute = int(obs.split()[7][2:4])
               second = float(obs.split()[8])
               pick_time = UTCDateTime(year, month, day, hour, minute, second)
               if obs.split()[4] == 'P':
                  ax.vlines(pick_time - streams[station][0].stats.starttime - time_diff, min(streams[station][0].data / np.amax(streams[station][0].data)) + distance, max(streams[station][0].data / np.amax(streams[station][0].data)) + distance, linestyles='solid', colors ='c', label='P-pick')
               if obs.split()[4] == 'S':
                  ax.vlines(pick_time - streams[station][0].stats.starttime - time_diff, min(streams[station][0].data / np.amax(streams[station][0].data)) + distance, max(streams[station][0].data / np.amax(streams[station][0].data)) + distance, linestyles='solid', colors ='r', label='S-pick')

    plt.xlabel('Seconds after event')
    plt.ylabel('Epicentral distance [km]')
    plt.grid(True)

    savename = glob.glob('*_tiebenn_loc')[0] + '/picks_location.pdf'

    plt.savefig(savename)
    plt.close()


def plot_hypoc_confidence_ellipsoid():
    """

    Generates a plot of the events location (maximum probability point), the probability density function (scatterplot) and the confidence ellipsoid as retrieved from the location with NonLinLoc.

    Args:
         None, because it works using files already produced during the location with NonLinLoc

    Returns:
         NLL_confidence_ellipsoid.pdf: the confidence ellipsoid projected on the XY, XZ and YZ planes and depicting the PDF scatterplot, as well as a red star indicating the point of maximum probability (hypocenter) and the 3 68% confidence ellipsoids projected on each plane
    """
    ls_file = glob.glob('gmt_*.LS.*')[0]
    e101_file = glob.glob('gmt_*LE_101*')[0]

    X = []; Y = []

    with open(glob.glob('*_tiebenn_loc/*scat.XY')[0], 'r') as f:
         for line in f:
             X.append(float(line.split()[0]))
             Y.append(float(line.split()[1]))

    ls_lines = open(ls_file, 'r').readlines()

    stations = []
    for i, line in enumerate(ls_lines):
        if line.startswith('# Station'):
           stations.append([ls_lines[i+2].split()[0], ls_lines[i+2].split()[1], ls_lines[i+2].split()[6]])

    with open(glob.glob('*_tiebenn_loc/loc*.hyp')[0], 'r') as f0:
         for line in f0:
             if len(line.split()) > 0:
                if line.split()[0] == 'HYPOCENTER':
                   max_x = float(line.split()[2])
                   max_y = float(line.split()[4])
                   max_z = float(line.split()[6])
    
    e101_lines = open(e101_file, 'r').readlines()

    ellipsoid1 = {}; ellipsoid2 = {}; ellipsoid3 = {}
    switch = 0

    region_xy = [max_x - 10, max_x + 10, max_y - 10, max_y + 10]

    if max_z <= 30.:
       region_zy = [0, 30, max_y - 10, max_y + 10]
       region_xz = [max_x - 10, max_x + 10, 0, 30]
    else:
         region_zy = [0, max_z, max_y - 10, max_y + 10]
         region_xz = [max_x - 10, max_x + 10, 0, max_z]

    for i, line in enumerate(e101_lines):
        if line.startswith('# Error'):
           count = 0; j = i + 3; e11 = []; e12 = []; e21 = []; e22 = []; e31 = []; e32 = []

           while count < 3:

                 if len(e101_lines[j].split()) == 2:
                    if count == 0:
                       e11.append(float(e101_lines[j].split()[0]))
                       e12.append(float(e101_lines[j].split()[1]))
                    if count == 1:
                       e21.append(float(e101_lines[j].split()[0]))
                       e22.append(float(e101_lines[j].split()[1]))
                    if count == 2:
                       e31.append(float(e101_lines[j].split()[0]))
                       e32.append(float(e101_lines[j].split()[1]))

                 if len(e101_lines[j].split()) == 1:
                    if e101_lines[j].split()[0] == '>':
                       count = count + 1

                 j = j + 1

           if switch == 0:
              ellipsoid1['e1'] = e11, e12
              ellipsoid1['e2'] = e21, e22
              ellipsoid1['e3'] = e31, e32
              switch = 1
           elif switch == 1:
              ellipsoid2['e1'] = e11, e12
              ellipsoid2['e2'] = e21, e22
              ellipsoid2['e3'] = e31, e32
              switch = 2
           else:
              ellipsoid3['e1'] = e11, e12
              ellipsoid3['e2'] = e21, e22
              ellipsoid3['e3'] = e31, e32

    pygmt.config(FONT='14p', FONT_ANNOT='12p')
    fig = pygmt.Figure()
    fig.basemap(region=region_xy, projection='X10c', frame=['xafg+lEast (km)', 'yafg+lNorth (km)', 'Wsen+tLocation: confidence ellipsoid'])
    fig.plot(x=X, y=Y, style='c1.25p', fill='black')

    for i in ellipsoid2:
        fig.plot(x=ellipsoid2[i][0], y=ellipsoid2[i][1], pen='1p,cyan')

    fig.plot(x=max_x, y=max_y, style='a6p', fill='red')

    for i in range(len(stations)):
        fig.text(x=stations[i][0], y=stations[i][1], text=stations[i][2])

    fig.shift_origin(xshift='w+0.5c')

    Z = []; Y = []
    with open(glob.glob('*_tiebenn_loc/*scat.ZY')[0], 'r') as f:
         for line in f:
             Z.append(float(line.split()[0]))
             Y.append(float(line.split()[1]))

    fig.basemap(region=region_zy, projection='X5c/10c', frame=['xafg+lDepth (km)', 'yafg+lNorth (km)', 'wSEn'])
    fig.plot(x=Z, y=Y, style='c1.25p', fill='black')

    for i in ellipsoid3:
        fig.plot(x=ellipsoid3[i][0], y=ellipsoid3[i][1], pen='1p,cyan')

    fig.plot(x=max_z, y=max_y, style='a6p', fill='red')
    fig.shift_origin(xshift='w-15.5c', yshift='h-15.5c')

    fig.basemap(region=region_xz, projection='X10c/-5c', frame=['xafg+lEast (km)', 'yafg+lDepth (km)', 'WSen'])

    X = []; Z = []
    with open(glob.glob('*_tiebenn_loc/*scat.XZ')[0], 'r') as f:
         for line in f:
             X.append(float(line.split()[0]))
             Z.append(float(line.split()[1]))

    fig.plot(x=X, y=Z, style='c1.25p', fill='black')

    for i in ellipsoid1:
        fig.plot(x=ellipsoid1[i][0], y=ellipsoid1[i][1], pen='1p,cyan')

    fig.plot(x=max_x, y=max_z, style='a6p', fill='red')

    with open(glob.glob('*_tiebenn_loc/loc_eqdatetime*.hyp')[0], 'r') as f:
         for line in f:
             if len(line.split()) > 0:
                lin_ = line.split()
                if lin_[0] == 'GEOGRAPHIC':
                   datetime = lin_[2] + '-' + lin_[3] + '-' + lin_[4] + ' ' + lin_[5] + ':' + lin_[6] + ':' + lin_[7]
                   geo_epi = 'Lat: ' + lin_[9] + '  Lon: ' + lin_[11]
                   depth = 'Depth: ' + lin_[13] + ' km'
                if lin_[0] == 'QUALITY':
                   rms_nphs = 'RMS: ' + lin_[8] + 's ' + lin_[9] + ':' + lin_[10]
                   gap_dist = 'Az_gap: ' + lin_[12] + ' Ne_sta: ' + lin_[14] + ' km'
                if lin_[0] == 'QML_ConfidenceEllipsoid':
                   ell_axes = 'semiMaj: ' + lin_[2] + ' semiInt: ' + lin_[6]
                   semi_min = 'semiMin: ' + lin_[4] 

    fig.text(x=11, y=10, text=datetime, font='11p,Courier-Bold,black', justify='LM', no_clip=True)
    fig.text(x=[11, 11], y=[12, 14], text=[geo_epi, depth], font='10p,Courier,black', justify='LT', no_clip=True)
    fig.text(x=[11, 11], y=[16, 18], text=[rms_nphs, gap_dist], font='10p,Courier,black', justify='LT', no_clip=True)
    fig.text(x=[11, 11, 11], y=[20, 22, 24], text=['NLL Confidence ellipsoid axes (km):', ell_axes, semi_min], font='10p,Courier,black', justify='LT', no_clip=True)

    fig.savefig(glob.glob('*_tiebenn_loc/')[0] + 'NLL_confidence_ellipsoid.pdf')

    for dele in glob.glob('gmt*.gmt'):
        os.remove(dele)

    loc_file = glob.glob('*_tiebenn_loc/loc_eqdatetime*.hyp')[0]
    new_name = glob.glob('*_tiebenn_loc')[0] + '/event_location.NLL'
    os.rename(loc_file, new_name)

    for dele in glob.glob('*_tiebenn_loc/loc_eqdatetime*'):
        os.remove(dele)

    return
