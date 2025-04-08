from importlib.resources import files
from math import floor

import numpy as np


def load_velmod(name):
    """

    Accesses the 1D velocity models in the TieBeNN package for seismic location with NonLinLoc.

    Args:
         model (int): The P- and S-wave velocity models to be loaded. The velocity models are read from files in the directory data/velocity_models/. New velocity models can be added following the NLL structure. See velocity models currently available in velmods
    Returns:
         velmod_list (list): Layers of a 1D velocity model as read from TieBeNN's resources
    """
    velmod = f"v{str(name)}"

    try:
        model_path = files('tiebenn.data.velocity_models').joinpath(velmod)

        velmod_list = []
        with model_path.open('r') as f:
             for line in f:
                 velmod_list.append(line.strip())
    except:
           class VelModLoadError(Exception):
                 pass
           raise VelModLoadError('Velocity model selected does not exist.')

    return velmod_list


def velmods(model, ev_lon, ev_lat):
    """
    
    Produces a seismic velocity model to be included in the NonLinLoc control file for hypocentre location. When no density model is associated to each model, density is calculated as a function of the P-wave velocity after Gardner et al. (1974) in g/cm^3. The density is, however, a carry over in the layer format from its original use in a waveform modelling code. It is not used in the NLL programs, so any convenient numerical value can be used.
    IMPORTANT: This functions is programmed so that each column will be separated by 3 blank spaces from the next one. This becomes of relevance during the phase association and must not be modified

    Args:
         model (int): The P- and S-wave velocity models to be exported. The velocity models are read from files in the directory data/velocity_models/. New velocity models can be added following the NLL structure. Currently, the following models are available:
           0 = IASP91
           1 = AK135 (for continental structure)
           2 = BGR velocity model (by J. Schlittenhardt)
           3 = WET, velocity model for Germany and adjacent regions
           4 = WB2012, Vogtland/ West Bohemia velocity model (1D, P+S), averaged from the 3D velocity models of Ruzek & Horalek (2013)
           5 = DEU, 3-layer model for Germany with a Moho depth of 28.5 km. It seems to be more appropriate for southern Germany. A vp/vs ratio of 1.68 was used to generate S-wave velocities
           6 = Crust1.0, for each latitude and longitude pair it produces up to 9 layers. Crust1.0 contains the topography at the epicenter
           7 = Crust1.0 + AK135, the Crust1.0 model is used for the crustal structure and the lower continental structure from AK135 is merged immediately below
           8 = Crustal velocity model from Küperkoch et al. (2018) for Insheim, (1D, P+S)
           9 = Regional model for LI (J. Borns)
          10 = Landau model (1D, P+S) for the northern Upper Rhine Graben
          11 = 1D model for Central Alps (Diehl et al. 2021)
          12 = WARNING! I STRONGLY SUSPECT THAT I HAVE NOT INTEGRATED THIS MODEL IN TIEBENN PROPERLY AND THE AUTHOR (DIEHL, ETH ZÜRICH) SHOULD BE CONTACTED FOR FURTHER CLARIFICATION!!(P+S) 3D model for Central Alps (Diehl et al. 2021)
          13 = 3D WEG model: high-resolution P-velocity model for the area around Rotenburg, Soehlingen, Soltau, Verden... (This model is not present in the current version of Tiebenn because it is really large and its implementation must be carefully tested!)
          14 = AlpsLocPS (Brazsus et al., 2024). 1D velocity models for the Greater Alpine Region
          15 = BENS (Reamer & Hinzen 2004): 1D P- and S-velocity models of the Northern Rhine Area
          16 = ASZmod1 (Mader et al. 2021): 1D P- and S-velocity models of the Albstadt Shear Zone
          17 = 3D P+S velocity models for the Upper Rhine Graben after Lengline et al. (2023)
          18 = PO1 (Malek et al. 2023): 1D P- and S-velocity model for the Novy-Kostel earthquake swarm region, constrained models
          19 = PO1 (Malek et al. 2023) + WB2012 from 11 km depth
         ev_lon (float): Event longitude. This parameter will be used if model = 6 or 7
         ev_lat (float): Event latitude. This parameter will be used if model = 6 or 7

    Returns:
         velmod (list): List of 1D velocity/density models in NonLinLoc format. The elements of the list are strings containing the information of each layer: LAYER, depth [km], p-wave velocity and velocity gradient, s-wave velocity and velocity gradient [km/s], density and density gradient in any convenient numerical value (see above)
    """
    try:
        if model not in [6, 7, 12, 13, 17]:
           velmod = load_velmod(model)

        if model == 6:
           crust1 = crustModel()
           velmod = []
           for i in crust1.get_point(ev_lat, ev_lon):
               velmod.append(f"LAYER   {str(-crust1.get_point(ev_lat, ev_lon)[i][4])}   {str(crust1.get_point(ev_lat, ev_lon)[i][0])}   0.0   {str(crust1.get_point(ev_lat, ev_lon)[i][1])}   0.0   {str(crust1.get_point(ev_lat, ev_lon)[i][2])}    0.0")

        if model == 7:
           crust1 = crustModel()
           velmod = []
           for i in crust1.get_point(ev_lat, ev_lon):
               velmod.append(f"LAYER   {str(-crust1.get_point(ev_lat, ev_lon)[i][4])}   {str(crust1.get_point(ev_lat, ev_lon)[i][0])}   0.0   {str(crust1.get_point(ev_lat, ev_lon)[i][1])}   0.0   {str(crust1.get_point(ev_lat, ev_lon)[i][2])}    0.0")
           velmod.append('LAYER   35.0   8.04   0.0   4.48   0.0   3.38   0.0')
           velmod.append('LAYER   77.5   8.045   0.0   4.49   0.0   3.38   0.0')
           velmod.append('LAYER   120.0   8.05   0.0   4.5   0.0   3.36   0.0')

    except:
           class WrongVelocityModelError(Exception):
                 pass
                 
           raise WrongVelocityModelError('Incorrect/undefined velocity model!')

    return velmod


class crustModel:
    """
    
    Top level model object to retreive information from the LLNL Crust 1.0 model. Created by J. Leeman and C. Ammon (https://github.com/jrleeman/Crust1.0).

    Args:
         vp (ndarray) P-wave velocity model
         vs (ndarray) S-wave velocity model
         rho (ndarray) Density model
         bnds (ndarray) Elevation of the top of the given layer with respect to sea level model
         layer_names (list) Names of the nine possible layers in the model
    """
    crust1_dir = files('tiebenn.data.crust1')

    def __init__(self):
        self.vp = np.loadtxt(crust1_dir.joinpath('crust1.vp'))
        self.vs = np.loadtxt(crust1_dir.joinpath('crust1.vs'))
        self.rho = np.loadtxt(crust1_dir.joinpath('crust1.rho'))
        self.bnds = np.loadtxt(crust1_dir.joinpath('crust1.bnds'))

        self.vp = self.vp.reshape((180, 360, 9))
        self.vs = self.vs.reshape((180, 360, 9))
        self.rho = self.rho.reshape((180, 360, 9))
        self.bnds = self.bnds.reshape((180, 360, 9))

        self.layer_names = ['water', 'ice', 'upper_sediments', 'middle_sediments', 'lower_sediments', 'upper_crust', 'middle_crust', 'lower_crust', 'mantle']

    def _get_index(self, lat, lon):
        """
        
        Returns in index values used to query the model for a given lat lon.

        Args:
             lat (float): Latitude of interest
             lat (float): Longitude of interest

        Returns:
             ilat (int) Index for given latitude
             ilon (int) Index for given longitude
        """

        if lon > 180:
            lon -= 360
        if lon < -180:
            lon += 360

        ilat = floor(90. - lat)
        ilon = floor(180 + lon)

        return int(ilat), int(ilon)

    def get_point(self, lat, lon):
        """
        
        Returns a model for a given latitude and longitude. Note that the model is only defined on a 1 degree grid starting at 89.5 and -179.5.

        Args:
        lat (float) Latitude of interest
        lon (float) Longitude of interest

        Returns:
        model_layers (dict): Dictionary of layers with the keys as layer names and the values as a list of vp, vs, density, layer thickness, and the top of the layer with respect to sea level
        """

        ilat, ilon = self._get_index(lat, lon)

        thickness = np.abs(np.ediff1d(self.bnds[ilat, ilon], to_end=[0]))

        model_layers = dict()

        for i, layer in enumerate(self.layer_names):
            vp = self.vp[ilat, ilon][i]
            vs = self.vs[ilat, ilon][i]
            rho = self.rho[ilat, ilon][i]
            bnd = self.bnds[ilat, ilon][i]
            layer_thickness = thickness[i]

            if layer_thickness >= 0.01 or layer == 'mantle':
                model_layers[layer] = [vp, vs, rho, layer_thickness, bnd]

        return model_layers


def select_velmod(ev_lat, ev_lon):
    """

    In case the velocity model was set to automatic, this function checks if there is a velocity model for the epicenter region which has demonstrated to produce better hypocenter estimations than the standard 2-layer BGR model. If not, the velocity model will be set as the BGR velocity model for event location.

    Args:
         ev_lat (float): Epicenter latitude
         ev_lon (float): Epicenter longitude
    Returns:
            velmod (int): The velocity model which will be used for event location within NonLinLoc.
    """
    if (ev_lon > 9.7 and ev_lon <= 13 and ev_lat >= 47 and ev_lat <= 48) or (ev_lon > 13 and ev_lon <= 15 and ev_lat >= 47 and ev_lat <= 49):
       velmod = 14
    elif ev_lon >= 6 and ev_lon <= 8 and ev_lat >= 50 and ev_lat <= 51:
         velmod = 15
    elif ev_lon >= 8 and ev_lon <= 10 and ev_lat > 48 and ev_lat < 49:
         velmod = 16
    elif ev_lon >= 6 and ev_lon <= 9.7 and ev_lat >= 47 and ev_lat <= 48:
         velmod = 11
    else:
         velmod = 2

    print('Selected velocity model for the input latitude and longitude: %s' % str(velmod))

    return velmod
