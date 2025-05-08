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
    Produces a seismic velocity model to be included in the NonLinLoc control file for hypocenter location.

    When no density model is associated to each model, density can be calculated as a function of the P-wave velocity after Gardner et al. (1974) in g/cm^3. The density is, however, a carryover in the layer format from its original use in a waveform modelling code. It is not used in the NLL programs, so any convenient numerical value can be used.

    Parameters
    ----------
    model : int
        The P- and S-wave velocity models to be exported. Velocity models are read from files in the directory `data/velocity_models/`. New models can be added following the NLL structure.

        The following model numbers are currently supported::

            0  = IASP91
            1  = AK135 (for continental structure)
            2  = BGR velocity model (by J. Schlittenhardt)
            3  = WET: Germany and adjacent regions
            4  = WB2012: Vogtland/West Bohemia 1D (P+S)
            5  = DEU: 3-layer Germany model (Moho at 28.5 km, vp/vs = 1.68)
            6  = Crust1.0: up to 9 layers per lat/lon point
            7  = Crust1.0 + AK135 hybrid model
            8  = Küperkoch et al. (2018) for Insheim (1D, P+S)
            9  = Regional model for LI (J. Borns)
           10  = Landau model (1D, P+S) for northern Upper Rhine Graben
           11  = Central Alps 1D (Diehl et al. 2021)
           12  = Not fully integrated and tested! 3D Central Alps (Diehl, ETH Zürich)
           13  = 3D WEG: high-res P model (not functional)
           14  = AlpsLocPS (Brazsus et al., 2024) — Greater Alpine Region
           15  = BENS (Reamer & Hinzen 2004): Northern Rhine Area
           16  = ASZmod1 (Mader et al. 2021): Albstadt Shear Zone
           17  = 3D URG (Lengline et al. 2023): Upper Rhine Graben (currently not functional)
           18  = (Malek et al. 2023): Novy-Kostel swarm
           19  = PO1 + WB2012 from 11 km depth
           20  = KIT6 (Ritter et al. 2024): 1D, P+S for East Eifel Volcanic Field

    ev_lon : float
        Event longitude. Used when model = 6 or 7.

    ev_lat : float
        Event latitude. Used when model = 6 or 7.

    Returns
    -------
    velmod : list of str
        List of 1D velocity/density layers in NonLinLoc format. Each element is a string representing one layer:
        LAYER depth[km] vp vp_grad vs vs_grad density density_grad

        Note: density values are placeholders and not used in NLL itself.
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
    def __init__(self):
        self.vp = np.loadtxt(files('tiebenn.data.crust1').joinpath('crust1.vp'))
        self.vs = np.loadtxt(files('tiebenn.data.crust1').joinpath('crust1.vs'))
        self.rho = np.loadtxt(files('tiebenn.data.crust1').joinpath('crust1.rho'))
        self.bnds = np.loadtxt(files('tiebenn.data.crust1').joinpath('crust1.bnds'))

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
    elif ev_lon >= 6.75 and ev_lon <= 8.2 and ev_lat >= 50 and ev_lat <= 50.75:
         velmod = 20
    else:
         velmod = 2

    print('Selected velocity model for the input latitude and longitude: %s' % str(velmod))

    return velmod
