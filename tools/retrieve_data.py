from obspy import UTCDateTime
from obspy.clients.fdsn import Client as client_fdsn
from obspy.clients.filesystem.sds import Client as client_sds
import os, glob, json, shutil

def make_station_list(stations_json, client_list, ev_lon, ev_lat, start_time, end_time, channel_list=[], filter_network=[], filter_station=[]):
    """
    
    Uses fdsn to find available stations in a specific geographical location and time period. This is a modified version of makeStationList of EQTransformer and it does not stop with an error if one of the selected client names does not have data available in the requested region and time.  

    Args:
         stations_json (str): Path of the json file that will be returned
         client_list (list): List of client names e.g. ['IRIS', 'SCEDC', 'USGGS'] from where to extract stations
         ev_lon (float): Longitude of the located event
         ev_lat (float): Latitude of the located event
         start_time (str): Start datetime for the beginning of the period. Format: yyyy-mm-dd hh:mm:ss.ss
         end_time (str): End datetime for the beginning of the period. Format: yyyy-mm-dd hh:mm:ss.ss
         channel_list (list, default=[]): A list containing the desired channel codes. Downloads will be limited to these channels based on priority. Defaults to [] --> all channels
         filter_network (list, default=[]): A list containing the network codes that need to be avoided
         filter_station (list, default=[]): A list containing the station names that need to be avoided

    Returns:
         stations_json: A dictionary containing relevant information for the available stations, such as geographical coordinates, elevation, available channels
    """
    
    min_lat = ev_lat - 1.7
    max_lat = ev_lat + 1.7
    min_lon = ev_lon - 1.7
    max_lon = ev_lon + 1.7

    station_list = {}
    for cl in client_list:
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

                          if len(channel_list) > 0:
                             chan_priority = [ch[:2] for ch in channel_list]

                             for chnn in chan_priority:
                                 if chnn in [ch[:2] for ch in new_chan]:
                                    new_chan = [ch for ch in new_chan if ch[:2] == chnn]

                          if len(new_chan) > 0 and (station not in station_list):
                             station_list[str(station)] = {'network': net, 'channels': list(set(new_chan)), 'coords': [lat, lon, elv], 'client': cl}

        except:
               pass

    json_dir = os.path.dirname(stations_json)
    if not os.path.exists(json_dir):
       os.makedirs(json_dir)
    with open(stations_json, 'w') as fp:
         json.dump(station_list, fp)

    return


def generate_station_database():
    """
    
    Python script for the creation of a dictionary storing station coordinates in a json file.
    """

    client_list = ['http://192.168.11.220:8080', 'LMU', 'GFZ', 'ETH', 'INGV', 'IPGP', 'NIEP', 'ODC', 'RESIF', 'RASPISHAKE', 'ORFEUS']

    stations_json = 'database_clients_stations.json'
    min_lat = 20
    max_lat = 60
    min_lon = -4
    max_lon = 25
    filter_network = []
    filter_station = []
    channel_list=['HHZ', 'HHE', 'HHN', 'EHZ', 'EHN', 'EHE', 'SHZ', 'SHN', 'SHE', 'HGZ', 'HGN', 'HGE', 'DNZ', 'DNE', 'DNN', 'HH1', 'HH2', 'EH1', 'EH2', 'SH1', 'SH2', 'BHZ', 'BHN', 'BHE', 'BH1', 'BH2']

    station_list = {}

    for cl in client_list:
         try:
            print('Creating station list for', cl)
            station_list[cl] = {}
            inventory = client_fdsn(cl).get_stations(minlatitude=min_lat, maxlatitude=max_lat, minlongitude=min_lon, maxlongitude=max_lon, starttime=UTCDateTime(), endafter=UTCDateTime(), level='channel')

            for ev in inventory:
                net = ev.code
                if net not in filter_network:
                   for st in ev:
                       station = st.code

                       if station not in filter_station:

                          elv = st.elevation
                          lat = st.latitude
                          lon = st.longitude
                          new_chan = [ch.code for ch in st.channels if ch.sample_rate >= 20.]

                          if len(channel_list) > 0:
                             chan_priority = [ch[:2] for ch in channel_list]
                             for chnn in chan_priority:
                                 if chnn in [ch[:2] for ch in new_chan]:
                                    new_chan = [ch for ch in new_chan if ch[:2] == chnn]

                          if len(new_chan) > 0 and (station not in station_list):
                             station_list[cl][str(station)] = {'network': net, 'channels': list(set(new_chan)), 'coords': [lat, lon, elv], 'client': cl}
         except:
               print('WARNING: Something went wrong with client', cl)
               pass

    station_list['BGR'] = station_list.pop('http://192.168.11.220:8080')
    for i in station_list['BGR']:
        station_list['BGR'][i]['client'] = 'BGR'

    with open(stations_json, 'w') as fp:
         json.dump(station_list, fp)
         
#    shutil.move(glob.glob(stations_json)[0], 'utils')

    return
