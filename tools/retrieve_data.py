import os

from joblib import Parallel, delayed
from obspy import UTCDateTime
from obspy.clients.fdsn import Client as client_fdsn
from obspy.clients.filesystem.sds import Client as client_sds


def process_client(cl, ev_lon, ev_lat, start_time, end_time, channel_list=[], filter_network=[], filter_station=[]):
    """

    Process a single client: fetch station inventory and filter based on criteria. Returns a dictionary of station info.

    Args:
         cl (str): Client name e.g. ['IRIS', 'SCEDC', 'USGGS'] from where to extract stations
         ev_lon (float): Longitude of the located event
         ev_lat (float): Latitude of the located event
         start_time (str): Start datetime for the beginning of the period. Format: yyyy-mm-dd hh:mm:ss.ss
         end_time (str): End datetime for the beginning of the period. Format: yyyy-mm-dd hh:mm:ss.ss
         channel_list (list, default=[]): A list containing the desired channel codes. Downloads will be limited to these channels based on priority. Defaults to [] --> all channels
         filter_network (list, default=[]): A list containing the network codes that need to be avoided
         filter_station (list, default=[]): A list containing the station names that need to be avoided

    Returns:
         client_station_list: A dictionary of the stations, if necessary
    """
    min_lat = ev_lat - 1.7
    max_lat = ev_lat + 1.7
    min_lon = ev_lon - 1.7
    max_lon = ev_lon + 1.7

    client_station_list = {}
    try:
        inventory = client_fdsn(cl).get_stations(minlatitude=min_lat, maxlatitude=max_lat, minlongitude=min_lon, maxlongitude=max_lon, starttime=UTCDateTime(start_time), endtime=UTCDateTime(end_time), endafter=UTCDateTime(start_time), level='channel')

        for ev in inventory:
            net = ev.code
            if net not in filter_network:
                for st in ev:
                    station = st.code
                    if station not in filter_station:
                        elv = st.elevation
                        lat = st.latitude
                        lon = st.longitude
                        new_chan = [ch.code for ch in st.channels]
                        if channel_list:
                            chan_priority = [ch[:2] for ch in channel_list]
                            for chnn in chan_priority:
                                if chnn in [ch[:2] for ch in new_chan]:
                                    new_chan = [ch for ch in new_chan if ch[:2] == chnn]
                        if new_chan and (station not in client_station_list):
                            client_station_list[str(station)] = {'network': net, 'channels': list(set(new_chan)), 'coords': [lat, lon, elv], 'client': cl}
    except Exception as e:
           pass

    return client_station_list


def make_station_list(client_list, ev_lon, ev_lat, start_time, end_time, channel_list=[], filter_network=[], filter_station=[]):
    """

    Uses fdsn to find available stations in a specific geographical location and time period.

    Args:
         client_list (list): List of client names e.g. ['IRIS', 'SCEDC', 'USGGS'] from where to extract stations
         ev_lon (float): Longitude of the located event
         ev_lat (float): Latitude of the located event
         start_time (str): Start datetime for the beginning of the period. Format: yyyy-mm-dd hh:mm:ss.ss
         end_time (str): End datetime for the beginning of the period. Format: yyyy-mm-dd hh:mm:ss.ss
         channel_list (list, default=[]): A list containing the desired channel codes. Downloads will be limited to these channels based on priority. Defaults to [] --> all channels
         filter_network (list, default=[]): A list containing the network codes that need to be avoided
         filter_station (list, default=[]): A list containing the station names that need to be avoided

    Returns:
         stations_list: A dictionary containing relevant information for the available stations, such as geographical coordinates, elevation, available channels
    """

    if len(client_list) <= int(os.cpu_count() / 2):
       njobs = len(client_list)
    else:
         njobs = int(os.cpu_count() / 2)

    results = Parallel(n_jobs=njobs, backend='threading')(
        delayed(process_client)(cl, ev_lon, ev_lat, start_time, end_time, channel_list, filter_network, filter_station)
        for cl in client_list)

    station_list = {}
    for res in results:
        station_list.update(res)

    return station_list
