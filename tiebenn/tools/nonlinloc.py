import glob
from importlib.resources import files
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd

from obspy import UTCDateTime

from nllgrid import NLLGrid

from .velocity_models import velmods
from .visualization import epic_sta_plot


def extract_3d_velmod(model_name: str, workdir: Path):
    """

    Creates a temporary directory with 3D velocity models copied from Tiebenn's resources.

    Args:
         model_name (str): Name of velocity model (e.g.: 'diehl', 'lengline', 'weg')
         workdir (pathlib.Path):
    Returns:
         dest (path): Path with the copied 3D velocity model files
    """
    model_dir = files('tiebenn.data.velocity_models').joinpath(model_name)

    dest = workdir / model_name
    dest.mkdir(exist_ok=True)

    for f in model_dir.iterdir():
        shutil.copy(f, dest)

    return dest


def inp_files_nlloc_sb(ev_lon, ev_lat, ev_time, data, nll3d, velmod, min_detections=3, verbosity=0):
    """

    Prepares the input files required by NonLinLoc for a successfully detected event with SeisBench.

    Args:
         ev_lon (float): Longitude of the located event with AD-Detektor
         ev_lat (float): Latitude of the located event with AD-Detektor
         ev_time (str): Origin time of the event. Format: yyyy-mm-dd hh:mm:ss.ss
         data (dict): A dictionary with information about the stations with picks
         nll3d (bool): If True, it will prepare the control files for a NonLinLoc depth estimation with a 3D velocity grid
         velmod (int): The seismic velocity model to be used in NonLinLoc for hypocentre location. See current options in velocity_models.py
         min_detections (int): It will define the minimum amount of stations on which the event was detected. If this minimum amount of detections is not achieved, the code will exit without calculating the hypocenter if the required minimum is not achieved. Value must be >= 3. Default = 3
         verbosity (int): Verbosity level of NonLinLoc. Possible options are -1 (silent), 0 (errors only), 1 (higher level warnings), >= 2 (low level warnings + information). Default is 0

    Returns:
         station_coordinates.txt: text file with information about all the stations with seismic detections following the GTSRCE/LATLON syntax
         obs_ttimes.obs: text file with the detected arrival times following the LOCFILES/NLLOC_OBS syntax
         nlloc_control.in: text file with the parameters for execution of NonLinLoc. A full description of each parameter can be found in https://github.com/alomax/NonLinLoc/blob/master/nlloc_sample/run/nlloc_sample.in
         nlloc_control_s.in: text file with the parameters for execution of Grid2Time for S-waves
    """

    if not len(glob.glob('*_tiebenn_loc')):
       print('No P- and S-wave arrival times were detected for this event. Stop.')
       return

    print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
    print('---------------------------------------------------------------')
    print('Preparing input files for hypocentre location with NonLinLoc...')
    print('---------------------------------------------------------------')
    print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

    if nll3d:
       if velmod not in [12, 13, 17]:
          try:
              for rem in glob.glob('*_tiebenn_loc/nll_control*'):
                  os.remove(rem)
              os.remove(glob.glob('*_tiebenn_loc/obs_ttimes.obs')[0]); os.remove(glob.glob('*_tiebenn_loc/station_coordinates.txt')[0])
          except:
                 pass

    if velmod == 13:
       vpvs = 1.74

    if velmod == 17:
       stations16 = []; stations17 =[]

    with open(f"{glob.glob('*_tiebenn_loc/')[0]}station_coordinates.txt", 'w') as file1:
         for stations in range(len(glob.glob('*_tiebenn_loc/csv_picks/*.csv'))):
             station = glob.glob('*_tiebenn_loc/csv_picks/*.csv')[stations].split('/')[2].split('_')[0]
             exptxt = f"GTSRCE {station}   LATLON   {str(data[station]['coords'][0])}  {str(data[station]['coords'][1])}  0.0  {str(0.001 * data[station]['coords'][2])}\n"

             file1.write(exptxt)

             if velmod == 17:
                lat_ = data[station]['coords'][0]; lon_ = data[station]['coords'][1]
                if lat_ >= 48.66166 - 0.898 and lat_ <= 48.66166 + 0.898 and lon_ >= 7.77386 - 1.429 and lon_ <= 7.77386 + 1.429:
                   stations17.append(exptxt)
                else:
                     stations16.append(exptxt)
    file1.close()

    if velmod == 17:

       if len(stations17) != 0:
          with open(f"{glob.glob('*_tiebenn_loc/')[0]}station_coordinates17.txt", 'w') as file17:
               for st17 in stations17:
                   file17.write(st17)
          file17.close()

       if len(stations16) != 0:
          with open(f"{glob.glob('*_tiebenn_loc/')[0]}station_coordinates16.txt", 'w') as file16:
               for st16 in stations16:
                   file16.write(st16)
          file16.close()

    def prob2uncer(x):
        """
        Since SeisBench does not export pick-uncertainties (as of December 2024) and one does not simply calculate uncertainties, the pick-probabilities are used as proxy for the uncertainties through a simple equation.
        """
        uncert = 0.02 / x
        return uncert

    with open(f"{glob.glob('*_tiebenn_loc/')[0]}obs_ttimes.obs", 'w') as file2:
         for obs in glob.glob('*_tiebenn_loc/csv_picks/*.csv'):
             dataframe = pd.read_csv(obs)
             for pick in range(len(dataframe)):
                 year = str(UTCDateTime(dataframe.arrival_time[pick]).year)
                 month = '%02d' % UTCDateTime(dataframe.arrival_time[pick]).month
                 day = '%02d' % UTCDateTime(dataframe.arrival_time[pick]).day
                 hour = '%02d' % UTCDateTime(dataframe.arrival_time[pick]).hour
                 minute = '%02d' % UTCDateTime(dataframe.arrival_time[pick]).minute
                 second = UTCDateTime(dataframe.arrival_time[pick]).second + 1e-6 * UTCDateTime(dataframe.arrival_time[pick]).microsecond
                 uncertainty = '{:.2e}'.format(prob2uncer(dataframe.probability[pick]))
                 observation = f"{dataframe.station[pick]}   ?   HHZ   ?   {dataframe.phase[pick]}   ?   {year}{month}{day}   {hour}{minute}   {str(second)}   GAU   {str(uncertainty)}   -1.00e+00   -1.00e+00   -1.00e+00   1\n"

                 file2.write(observation)
    file2.close()

    with open(f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control.in", 'w') as file3:
         file3.write(f"CONTROL   {str(verbosity)}   54321\n")

         if velmod != 13:
            file3.write(f"TRANS   LAMBERT   Clarke-1880   {str(ev_lat)}   {str(ev_lon)}   {str(ev_lat - 1)}   {str(ev_lat + 1)}   {str(0.0)}\n")
         else:
              file3.write(f"TRANS   LAMBERT   WGS-84   52.3348   9.00   52    53    0.0\n")

         if not nll3d:
            if velmod not in [12, 13]:
               file3.write(f"VGOUT   {glob.glob('*_tiebenn_loc/')[0]}model_layer\n")
               file3.write(f"VGTYPE   P\n")
               file3.write(f"VGTYPE   S\n")
               file3.write(f"VGGRID   2   401   121   0.0   0.0   -3.0   1   1   1   SLOW_LEN\n")

               if velmod != 17:
                  for i in velmods(velmod, ev_lon, ev_lat):
                      file3.write(f"{i}\n")
               else:
                    for i in velmods(model=10, ev_lon=ev_lon, ev_lat=ev_lat):
                        file3.write(f"{i}\n")

            if velmod == 12:
               model_path = extract_3d_velmod(model_name='diehl', workdir='temp_velmodfiles')

               file3.write(f"VGINP    {model_path}/diehl_3D_pg.txt\n")
               file3.write(f"VGOUT   {glob.glob('*_tiebenn_loc/')[0]}model_layer\n")
               file3.write(f"VGTYPE   P\n")
               file3.write(f"VGGRID   43   43   13   -210.   -210   -7.0   10   10   10   SLOW_LEN   SLOW_LEN\n")

         if velmod not in [13, 17]:
            file3.write(f"GTFILES   {glob.glob('*_tiebenn_loc/')[0]}model_layer   {glob.glob('*_tiebenn_loc/')[0]}time_layer   P\n")
         elif velmod == 13:
              model_path = extract_3d_velmod(model_name='weg', workdir='temp_velmodfiles')

              file3.write(f"GTFILES   {model_path}/GitterWEG   {glob.glob('*_tiebenn_loc/')[0]}time_layer   P\n")

         if not nll3d:
            if velmod not in [12, 13, 17]:
               file3.write(f"GTMODE   GRID2D   ANGLES_NO\n")
            elif velmod in [12, 13]:
                 file3.write(f"GTMODE   GRID3D   ANGLES_NO\n")
         else:
              file3.write(f"GTMODE   GRID3D   ANGLES_NO\n")

         if velmod != 17:
            file3.write(f"INCLUDE   {glob.glob('*_tiebenn_loc/')[0]}station_coordinates.txt\n")
            file3.write(f"GT_PLFD   1.0e-3   0\n")

         file3.write(f"LOCSIG   Tiebenn/NonLinLoc\n")
         file3.write(f"LOCCOM   {str(ev_time)}\n")
         file3.write(f"LOCFILES   {glob.glob('*_tiebenn_loc/')[0]}obs_ttimes.obs   NLLOC_OBS   {glob.glob('*_tiebenn_loc/')[0]}time_layer   {glob.glob('*_tiebenn_loc/')[0]}loc_eqdatetime\n")
         file3.write(f"LOCHYPOUT   SAVE_NLLOC_ALL   NLL_FORMAT_VER_2\n")
         file3.write(f"LOCSEARCH   OCT   10   10   10   0.01   20000   9000   0   1\n")

         if velmod == 13:
            file3.write(f"LOCGRID   884   733   150    4.15   31.55  -0.2  0.1 0.1 0.1  PROB_DENSITY   SAVE\n")
            file3.write(f"LOCMETH   EDT_OT_WT   9999.0   3   -1   0   {str(vpvs)}   -1   0   1\n")
         elif velmod == 17:
              file3.write(f"LOCGRID   401   401   101   -100   -100   -1.0   0.5   0.5   0.5  PROB_DENSITY   SAVE\n")
              file3.write(f"LOCMETH   EDT_OT_WT   9999.0   3   -1   0   -1   -1   0   1\n")
         else:
              file3.write(f"LOCGRID   201   201   101   -100.0   -100.0   0.0   1   1   1   PROB_DENSITY   SAVE\n")
              file3.write(f"LOCMETH   EDT_OT_WT   9999.0   3   -1   0   -1   -1   0   1\n")

         file3.write(f"LOCGAU   0.2   0.0\n")
         file3.write(f"LOCGAU2   0.02   0.05   2.0\n")
         file3.write(f"LOCQUAL2ERR   0.1   0.5   1.0   2.0   99999.9\n")
         file3.write(f"LOCANGLES   ANGLES_NO   5\n")
         file3.write(f"LOCPHSTAT   9999.0   -1   9999.0   1.0   1.0   9999.9   -9999.9   9999.9\n")
    file3.close()

    if velmod != 17:
       with open(f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control_s.in", 'w') as file4:
            file4.write(f"CONTROL   {str(verbosity)}   54321\n")

            if velmod != 13:
               file4.write(f"TRANS   LAMBERT   Clarke-1880   {str(ev_lat)}   {str(ev_lon)}   {str(ev_lat - 1)}   {str(ev_lat + 1)}    {str(0.0)}\n")
            else:
                 file4.write(f"TRANS   LAMBERT   WGS-84   52.3348   9.00   52    53    0.0\n")

            if velmod == 12:
               model_path = extract_3d_velmod(model_name='diehl', workdir='temp_velmodfiles')

               file4.write(f"VGINP    {model_path}/diehl_3D_sg.txt\n")
               file4.write(f"VGOUT   {glob.glob('*_tiebenn_loc/')[0]}model_layer\n")
               file4.write(f"VGTYPE   S\n")
               file4.write(f"VGGRID   43   43   13   -210.   -210   -7.0   10   10   10   SLOW_LEN   SLOW_LEN\n")

            if velmod != 13:
               file4.write(f"GTFILES   {glob.glob('*_tiebenn_loc/')[0]}model_layer   {glob.glob('*_tiebenn_loc/')[0]}time_layer   S\n")
            else:
                 model_path = extract_3d_velmod(model_name='weg', workdir='temp_velmodfiles')

                 file4.write(f"GTFILES   {model_path}/GitterWEG   {glob.glob('*_tiebenn_loc/')[0]}time_layer   S\n")

            if not nll3d:
               if velmod not in [12, 13]:
                  file4.write(f"GTMODE   GRID2D   ANGLES_NO\n")
               else:
                    file4.write(f"GTMODE   GRID3D   ANGLES_NO\n")
            else:
                 file4.write(f"GTMODE   GRID3D   ANGLES_NO\n")

            file4.write(f"INCLUDE   {glob.glob('*_tiebenn_loc/')[0]}station_coordinates.txt\n")
            file4.write(f"GT_PLFD   1.0e-3   0\n")
       file4.close()

    if velmod == 17:
       if len(stations17) != 0:
          model_path = extract_3d_velmod(model_name='lengline', workdir='temp_velmodfiles')

          with open(f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control17.in", 'w') as file317:
               file317.write(f"CONTROL   {str(verbosity)}   54321\n")
               file317.write(f"TRANS    LAMBERT WGS-84   48.66166   7.77386   47   49   0.0\n")
               file317.write(f"GTFILES   {model_path}/layer   {glob.glob('*_tiebenn_loc/')[0]}time_layer   P\n")
               file317.write(f"GTMODE   GRID3D   ANGLES_NO\n")
               file317.write(f"INCLUDE   {glob.glob('*_tiebenn_loc/')[0]}station_coordinates17.txt\n")
               file317.write(f"GT_PLFD   1.0e-3   0\n")
          file317.close()

          with open(f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control17_s.in", 'w') as file417:
               file417.write(f"CONTROL   {str(verbosity)}   54321\n")
               file417.write(f"TRANS    LAMBERT WGS-84   48.66166   7.77386   47   49   0.0\n")
               file417.write(f"GTFILES   {model_path}/layer   {glob.glob('*_tiebenn_loc/')[0]}time_layer   S\n")
               file417.write(f"GTMODE   GRID3D   ANGLES_NO\n")
               file417.write(f"INCLUDE   {glob.glob('*_tiebenn_loc/')[0]}station_coordinates17.txt\n")
               file417.write(f"GT_PLFD   1.0e-3   0\n")
          file417.close()

       if len(stations16) != 0:
          with open(f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control16.in", 'w') as file316:
               file316.write(f"CONTROL   {str(verbosity)}   54321\n")
               file316.write(f"TRANS   LAMBERT   Clarke-1880   {str(ev_lat)}   {str(ev_lon)}   {str(ev_lat - 1)}   {str(ev_lat + 1)}   {str(0.0)}\n")
               file316.write(f"GTFILES    {glob.glob('*_tiebenn_loc/')[0]}model_layer    {glob.glob('*_tiebenn_loc/')[0]}time_layer   P\n")
               file316.write(f"GTMODE   GRID2D   ANGLES_NO\n")
               file316.write(f"INCLUDE   {glob.glob('*_tiebenn_loc/')[0]}station_coordinates16.txt\n")
               file316.write(f"GT_PLFD   1.0e-3   0\n")
          file316.close()

          with open(f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control16_s.in", 'w') as file416:
               file416.write(f"CONTROL   {str(verbosity)}   54321\n")
               file416.write(f"TRANS   LAMBERT   Clarke-1880   {str(ev_lat)}   {str(ev_lon)}   {str(ev_lat - 1)}   {str(ev_lat + 1)}   {str(0.0)}\n")
               file416.write(f"GTFILES    {glob.glob('*_tiebenn_loc/')[0]}model_layer    {glob.glob('*_tiebenn_loc/')[0]}time_layer   S\n")
               file416.write(f"GTMODE   GRID2D   ANGLES_NO\n")
               file416.write(f"INCLUDE   {glob.glob('*_tiebenn_loc/')[0]}station_coordinates16.txt\n")
               file416.write(f"GT_PLFD   1.0e-3   0\n")
          file416.close()

    if os.path.exists('temp_velmodfiles'):
       shutil.rmtree('temp_velmodfiles')

    print('Files prepared.')


def pynlloc(control_file, control_file_s, data, velmod, nll3d, plots):
    """

    It runs the necessary executables for hypocenter estimation using NonLinLoc: Vel2Grid, Grid2Time and NLLoc.

    Args:
         control_file (str): the complete path and filename for the input control file with certain required and optional statements specifying program parameters and input/output file names and locations. These parameters and files are used to execute NonLinLoc
         control_file_s (str): the portion of the control file with the required statements, parameters and file locations for the execution of Grid2Time for S-waves. In the current version of NonLinLoc, the generation of traveltimes is under a non-repeatable parameter (GTFILES). As a consecuence, Grid2Time must be executed twice for generation of P- and S-wave traveltimes. Changing GTFILES to repeatable is as of September 2024 in the to-do list of the author A. Lomax (pers. commun.)
         velmod (int): The seismic velocity model to be used in NonLinLoc for hypocentre location. See current options in velmods.py
         data (dict): Dictionary with information of seismic stations
         nll3d (bool): If True, it will skip the execution of Vel2Grid
         plots (bool): if True, it will plot the epicenter with the stations used for location, if the location was sucessful

    Returns:
         event_location.NLL: here are the location and quality, statistical values in case the location was sucessful   
         epicenter_stations.pdf: map depicting the epicenter location on the map, alongside the stations with phase picks used by NonLinLoc for depth estimation
    """

    def Vel2Grid(control_file):
        """

        Subprocess to call and execute Vel2Grid, which converts analytic or other velocity model specifications into a 3D Grid file containing velocity or slowness values.

        Args:
             control_file (str): the control file to execute Vel2Grid

        Returns:
             Grid files containing velocity values. These files are deleted after pynlloc is executed
        """

        process = subprocess.Popen(['Vel2Grid', control_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        if len(err) > 0:
           print(err)

    def Vel2Grid3D(control_file):
        """

        Subprocess to call and execute Vel2Grid3D, which converts analytic or other 3D velocity model specifications into a 3D Grid file containing velocity or slowness values.

        Args:
             control_file (str): the control file to execute Vel2Grid3D

        Returns:
             Grid files containing velocity (slowness*lenght) values. These files are deleted after pynlloc is executed
        """

        process = subprocess.Popen(['Vel2Grid3D', control_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        if len(err) > 0:
           print(err)

    def Grid2Time(control_file):
        """

        Subprocess to call and execute Grid2Time, which, given a velocity model grid (previously obtained using Vel2Grid), calculates the travel-time from a source point in the grid to all other points in the grid.

        Args:
             control_file (str): the control file to execute Grid2Time

        Returns:
             Travel-times throughout a grid, which are written to a separate grid file for each phase at each station. These files are deleted after pynlloc is executed
        """

        process = subprocess.Popen(['Grid2Time', control_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        if len(err) > 0:
           print(err)

    def NLLoc(control_file):
        """

        Subprocess to call and execute NLLoc, which will produce an estimate of the posterior probability density function (PDF) for the spatial, x,y,z hypocenter location.

        Args:
             control_file (str): the control file to execute NLLoc

        Returns:
             Text files with location results. Most of these files are deleted after pynlloc is executed
        """

        process = subprocess.Popen(['NLLoc', control_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        if len(err) > 0:
           print(err)

    def Grid2GMT(control_file, gridloc, gmt_dir, plot_type1, plot_type2):
        """

        Subprocess to call and execute Grid2GMT, which converts NonLinLoc Grids to be read for location plots with PyGMT (GMT). More information and examples of using Grid2GMT in http://alomax.free.fr/nlloc/

        Args:
             control_file (str): the control file to execute Grid2GMT
             gridloc (str): the full path to the three grid files generated by NonLinLoc without extension (.hyp, .hdr, .scat)
             gmt_dir (str): full path to the directory where the exported GMT scripts will be saved
             plot_type1 (str): if velocity/slowness model, contour lines, error ellipsoids... In our case we use L mode (location)
             plot_type2 (str): two types of plots are relevant: S (scatterplot) and E101 (maximum likelihood point and confidence ellipsoid will be produced)
        Returns:
             - loc_eqdatetime.<date>.<time>.grid0.loc.scat.<XY/XZ/YZ>: The probability density function (PDF) scatterplot projected on the xy, xz and xz planes
             - <gmt_dir>loc_eqdatetime.<date>.<time>.grid0.loc.LE_101.gmt (if plot_type 2 = E101)
             - <gmt_dir>loc_eqdatetime.<date>.<time>.grid0.loc.LS.gmt (if plot_type 2 = S)
        """

        process = subprocess.Popen(['Grid2GMT', control_file, gridloc, gmt_dir, plot_type1, plot_type2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        if len(err) > 0:
           print(err)

    if not glob.glob(f"{glob.glob('*_tiebenn_loc/')[0]}nlloc_control.in"):
       print('NLL: Stop.')
       return

    print('NonLinLoc: Locating event...')

    try:
        if not nll3d:
           print('NLL: Producing velocity/slowness model from velocity description...')
           if velmod not in [12, 13]:
              Vel2Grid(control_file)
           elif velmod == 12:
              Vel2Grid3D(control_file)
              Vel2Grid3D(control_file_s)

        if not nll3d:
           print('NLL: Calculating P-wave travel times...')
        else:
             print('NLL3D: Calculating P-wave travel times...')

        if velmod != 17:
           Grid2Time(control_file)
        else:
             for control in glob.glob('*_tiebenn_loc/nlloc_control17*'):
                 Grid2Time(control)

             for rem_p in glob.glob(f"{glob.glob('*_tiebenn_loc/')[0]}time_layer.P.mod*"):
                 os.remove(rem_p)
             for rem_s in glob.glob(f"{glob.glob('*_tiebenn_loc/')[0]}time_layer.S.mod*"):
                 os.remove(rem_s)

        if not nll3d:
           print('NLL: Calculating S-wave travel times...')
        else:
             print('NLL3D: Calculating S-wave travel times...')

        if velmod != 17:
           Grid2Time(control_file_s)
        else:
             for control in glob.glob('*_tiebenn_loc/nlloc_control16*'):
                 Grid2Time(control)

        if not nll3d:
           print('NLL: Depth estimation...')
        else:
             print('NLL3D: Depth estimation...')
        NLLoc(control_file)

        location = 'Location completed.'
        location_check = False
        reject = 'REJECTED'
        with open(glob.glob('*_tiebenn_loc/*.hyp')[0], 'r') as f:
             lines = f.readlines()
             for line in lines:
                 if line.find(reject) != -1:
                    print('NonLinLoc: location was', reject)
                 else:
                      if line.find(location) != -1:
                         location_check = True
                         print('##############################################################')
                         if not nll3d:
                            print('NonLinLoc:', location)
                         else:
                              print('NonLinLoc3D:', location)
                         for out_line in lines:
                             if out_line.find('GEOGRAPHIC') != -1:
                                print('Origin time:', f"{out_line.split()[4]}-{out_line.split()[3]}-{out_line.split()[2]}_{out_line.split()[5]}:{out_line.split()[6]}:{out_line.split()[7]}")
                                print('            ', f"{out_line.split()[8]}:{out_line.split()[9]}", f"{out_line.split()[10]}:{out_line.split()[11]}", f"{out_line.split()[12]}:", out_line.split()[13], 'km')
                                event_latitude = float(out_line.split()[9])
                                event_longitude = float(out_line.split()[11])
                                print('--------------------------------------------------------------')
                                print('Location quality:')
                             if out_line.find('STATISTICS') != -1:
                                print('Ellipsoid semi-major axis:', out_line.split()[-1])
                             if out_line.find('QUALITY') != -1:
                                print(f"{out_line.split()[7]}:{out_line.split()[8]}", 'Number of phases:', out_line.split()[10], f"{out_line.split()[11]}:", out_line.split()[12], 'Distance from hypocenter to nearest station:', out_line.split()[14], 'km')
                                sta_gap = float(out_line.split()[12])

                         print('##############################################################')

             for rem1 in glob.glob('*_tiebenn_loc/last*'):
                 os.remove(rem1)
             for rem2 in glob.glob('*_tiebenn_loc/time_layer*'):
                 os.remove(rem2)
             for rem3 in glob.glob('*_tiebenn_loc/model_*'):
                 os.remove(rem3)
             for rem4 in glob.glob('*_tiebenn_loc/*.sum.*'):
                 os.remove(rem4)
             for rem5 in glob.glob('*_tiebenn_loc/*eqdatetime_nlloc*'):
                 os.remove(rem5)

    except:
          class NLLError(Exception):
                pass
          raise NLLError('NonLinLoc: Something went wrong. Event was not located')

    sta_nearest = 9999.
    if location_check:
       stations_file = f"{glob.glob('*_tiebenn_loc/')[0]}station_coordinates.txt"
       with open(stations_file) as f:
            for line in f:
                sta_dist_ = float(data[line.split()[1]]['epic_distance'])
                if sta_dist_ <= sta_nearest:
                   sta_nearest = sta_dist_

       if plots:
                try:
                    epic_sta_plot()
                except:
                       print('Something went wrong with the GMT plot. Was the installation of GMT and/or PyGMT correct?')
                       pass

                gridloc = glob.glob('*_tiebenn_loc/loc_eqdatetime.*.*.grid0.loc.hyp')[0][:-4]
                Grid2GMT(control_file, gridloc, 'gmt_', 'L', 'S')
                Grid2GMT(control_file, gridloc, 'gmt_', 'L', 'E101')
    else:
         sta_gap = 9999.; sta_nearest = 9999.

    return sta_gap, sta_nearest


def create3dgrid(ev_lon, ev_lat, velmod):
    """

    This function uses NLLGrid (https://github.com/claudiodsf/nllgrid) to transform 3D velocity descriptions (e.g. the Crust1 model, which is implemented in Tiebenn) into a 3D SLOW_LEN grid for the calculation of travel times with NonLinLoc. Assumption: 1° equals 111 km.

    Args:
         ev_lon (float): Longitude of the located event with AD-Detektor
         ev_lat (float): Latitude of the located event with AD-Detektor
         velmod (int): The seismic velocity model to be used in NonLinLoc for hypocentre location. See current options in velmods
    Returns:
         model_layer.(P/S).mod(.hdr/.buf): The 3D SLOW_LEN grids for P- and S-wave velocities
    """
    print('Creating 3D NLLGrids...')

    zmin = -1.0; zmax = 120.0; ymin = -277.5; ymax = 277.5; xmin = -277.5; xmax = 277.5
    dz = 1; dx = 1; dy = 1

    zvalues = np.linspace(start=zmin, stop=zmax, num=int(np.ceil((zmax - zmin)/dz)))

    # V33
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon, ev_lat):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v33 = {}; v33['vp'] = np.asarray(vp_values); v33['vs'] = np.asarray(vs_values)

    # V32
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 1, ev_lat):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v32 = {}; v32['vp'] = np.asarray(vp_values); v32['vs'] = np.asarray(vs_values)

    # V34
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 1, ev_lat):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v34 = {}; v34['vp'] = np.asarray(vp_values); v34['vs'] = np.asarray(vs_values)

    # V22
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 1, ev_lat + 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v22 = {}; v22['vp'] = np.asarray(vp_values); v22['vs'] = np.asarray(vs_values)

    # V23
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon, ev_lat + 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v23 = {}; v23['vp'] = np.asarray(vp_values); v23['vs'] = np.asarray(vs_values)

    # V24
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 1, ev_lat + 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v24 = {}; v24['vp'] = np.asarray(vp_values); v24['vs'] = np.asarray(vs_values)

    # V42
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 1, ev_lat - 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v42 = {}; v42['vp'] = np.asarray(vp_values); v42['vs'] = np.asarray(vs_values)

    # V43
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon, ev_lat - 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v43 = {}; v43['vp'] = np.asarray(vp_values); v43['vs'] = np.asarray(vs_values)

    # V44
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 1, ev_lat - 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v44 = {}; v44['vp'] = np.asarray(vp_values); v44['vs'] = np.asarray(vs_values)

    # V11
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 2, ev_lat + 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v11 = {}; v11['vp'] = np.asarray(vp_values); v11['vs'] = np.asarray(vs_values)

    # V12
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 1, ev_lat + 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v12 = {}; v12['vp'] = np.asarray(vp_values); v12['vs'] = np.asarray(vs_values)

    # V13
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon, ev_lat + 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v13 = {}; v13['vp'] = np.asarray(vp_values); v13['vs'] = np.asarray(vs_values)

    # V14
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 1, ev_lat + 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v14 = {}; v14['vp'] = np.asarray(vp_values); v14['vs'] = np.asarray(vs_values)

    # V15
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 2, ev_lat + 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v15 = {}; v15['vp'] = np.asarray(vp_values); v15['vs'] = np.asarray(vs_values)

    # V21
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 2, ev_lat + 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v21 = {}; v21['vp'] = np.asarray(vp_values); v21['vs'] = np.asarray(vs_values)

    # V25
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 2, ev_lat + 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v25 = {}; v25['vp'] = np.asarray(vp_values); v25['vs'] = np.asarray(vs_values)

    # V31
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 2, ev_lat):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v31 = {}; v31['vp'] = np.asarray(vp_values); v31['vs'] = np.asarray(vs_values)

    # V35
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 2, ev_lat):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v35 = {}; v35['vp'] = np.asarray(vp_values); v35['vs'] = np.asarray(vs_values)

    # V41
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 2, ev_lat - 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v41 = {}; v41['vp'] = np.asarray(vp_values); v41['vs'] = np.asarray(vs_values)

    # V45
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 2, ev_lat - 1):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v45 = {}; v45['vp'] = np.asarray(vp_values); v45['vs'] = np.asarray(vs_values)

    # V51
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 2, ev_lat - 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v51 = {}; v51['vp'] = np.asarray(vp_values); v51['vs'] = np.asarray(vs_values)

    # V52
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon - 1, ev_lat - 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v52 = {}; v52['vp'] = np.asarray(vp_values); v52['vs'] = np.asarray(vs_values)

    # V53
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon, ev_lat - 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v53 = {}; v53['vp'] = np.asarray(vp_values); v53['vs'] = np.asarray(vs_values)

    # V54
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 1, ev_lat - 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v54 = {}; v54['vp'] = np.asarray(vp_values); v54['vs'] = np.asarray(vs_values)

    # V55
    dep_list = []; vp_list_values = []; vs_list_values = []
    for layer in velmods(velmod, ev_lon + 2, ev_lat - 2):
        dep_list.append(float(layer.split()[1]))
        vp_list_values.append(float(layer.split()[2]))
        vs_list_values.append(float(layer.split()[4]))
    if float(layer.split()[1]) < 120.:
       dep_list.append(120.)
       vp_list_values.append(float(layer.split()[2]))
       vs_list_values.append(float(layer.split()[4]))

    index = 0; vp_values = []; vs_values = []
    for zstep in range(len(zvalues)):
        if zvalues[zstep] >= dep_list[index+1]:
           index = index + 1
        vp_values.append(vp_list_values[index])
        vs_values.append(vs_list_values[index])

    v55 = {}; v55['vp'] = np.asarray(vp_values); v55['vs'] = np.asarray(vs_values)

    yvalues = np.linspace(start=ymin, stop=ymax, num=int(np.ceil((ymax - ymin)/dy)))

    vp_y1 = []; vp_y2 = []; vp_y3 = []; vp_y4 = []; vp_y5 = []; vs_y1 = []; vs_y2 = []; vs_y3 = []; vs_y4 = []; vs_y5 = []

    for ystep in range(len(yvalues)):
        if yvalues[ystep] < -166.5:
           vp_y1.append(v51['vp']); vp_y2.append(v52['vp']); vp_y3.append(v53['vp']); vp_y4.append(v54['vp']); vp_y5.append(v55['vp'])
           vs_y1.append(v51['vs']); vs_y2.append(v52['vs']); vs_y3.append(v53['vs']); vs_y4.append(v54['vs']); vs_y5.append(v55['vs'])
        elif yvalues[ystep] >= -166.5 and yvalues[ystep] < -55.5:
           vp_y1.append(v41['vp']); vp_y2.append(v42['vp']); vp_y3.append(v43['vp']); vp_y4.append(v44['vp']); vp_y5.append(v45['vp'])
           vs_y1.append(v41['vs']); vs_y2.append(v42['vs']); vs_y3.append(v43['vs']); vs_y4.append(v44['vs']); vs_y5.append(v45['vs'])
        elif yvalues[ystep] >= -55.5 and yvalues[ystep] < 55.5:
           vp_y1.append(v31['vp']); vp_y2.append(v32['vp']); vp_y3.append(v33['vp']); vp_y4.append(v34['vp']); vp_y5.append(v35['vp'])
           vs_y1.append(v31['vs']); vs_y2.append(v32['vs']); vs_y3.append(v33['vs']); vs_y4.append(v34['vs']); vs_y5.append(v35['vs'])
        elif yvalues[ystep] >= 55.5 and yvalues[ystep] < 166.5:
           vp_y1.append(v21['vp']); vp_y2.append(v22['vp']); vp_y3.append(v23['vp']); vp_y4.append(v24['vp']); vp_y5.append(v25['vp'])
           vs_y1.append(v21['vs']); vs_y2.append(v22['vs']); vs_y3.append(v23['vs']); vs_y4.append(v24['vs']); vs_y5.append(v25['vs'])
        else:
           vp_y1.append(v11['vp']); vp_y2.append(v12['vp']); vp_y3.append(v13['vp']); vp_y4.append(v14['vp']); vp_y5.append(v15['vp'])
           vs_y1.append(v11['vs']); vs_y2.append(v12['vs']); vs_y3.append(v13['vs']); vs_y4.append(v14['vs']); vs_y5.append(v15['vs'])

    xvalues = np.linspace(start=xmin, stop=xmax, num=int(np.ceil((xmax - xmin)/dx)))

    vp_x = []; vs_x = []
    for xstep in range(len(xvalues)):
        if xvalues[xstep] < -166.5:
           vp_x.append(np.asarray(vp_y1))
           vs_x.append(np.asarray(vs_y1))
        elif xvalues[xstep] >= -166.5 and xvalues[xstep] < -55.5:
           vp_x.append(np.asarray(vp_y2))
           vs_x.append(np.asarray(vs_y2))
        elif xvalues[xstep] >= -55.5 and xvalues[xstep] < 55.5:
           vp_x.append(np.asarray(vp_y3))
           vs_x.append(np.asarray(vs_y3))
        elif xvalues[xstep] >= 55.5 and xvalues[xstep] < 166.5:
           vp_x.append(np.asarray(vp_y4))
           vs_x.append(np.asarray(vs_y4))
        else:
           vp_x.append(np.asarray(vp_y5))
           vs_x.append(np.asarray(vs_y5))

    vel_p = dx / np.asarray(vp_x)
    vel_s = dx / np.asarray(vs_x)

    grd_p = NLLGrid()
    grd_s = NLLGrid()

    grd_p.array = vel_p
    grd_p.dx = dx
    grd_p.dy = dy
    grd_p.dz = dz
    grd_p.x_orig = xmin
    grd_p.y_orig = ymin
    grd_p.z_orig = zmin
    grd_p.type = 'SLOW_LEN'
    grd_p.orig_lat = ev_lat
    grd_p.orig_lon = ev_lon
    grd_p.proj_name = 'LAMBERT'
    grd_p.first_std_paral = ev_lat - 1
    grd_p.second_std_paral = ev_lat + 1
    grd_p.proj_ellipsoid = 'Clarke-1880'
    grd_p.basename = 'model_layer.P.mod'
    grd_p.write_hdr_file()
    grd_p.write_buf_file()

    grd_s.array = vel_s
    grd_s.dx = dx
    grd_s.dy = dy
    grd_s.dz = dz
    grd_s.x_orig = xmin
    grd_s.y_orig = ymin
    grd_s.z_orig = zmin
    grd_s.type = 'SLOW_LEN'
    grd_s.orig_lat = ev_lat
    grd_s.orig_lon = ev_lon
    grd_s.proj_name = 'LAMBERT'
    grd_s.first_std_paral = ev_lat - 1
    grd_s.second_std_paral = ev_lat + 1
    grd_s.proj_ellipsoid = 'Clarke-1880'
    grd_s.basename = 'model_layer.S.mod'
    grd_s.write_hdr_file()
    grd_s.write_buf_file()

    try:
        for grid in glob.glob('model_layer*'):
            shutil.move(grid, glob.glob('*_tiebenn_loc/')[0])
    except:
           print('Problem moving grid to *_tiebenn_loc directory.')
           pass
    return
