# Velocity models implemented in TieBeNN

This is a list of the velocity models currently implemented in TieBeNN. As further local, dedicated velocity models are published/discovered, they can be integrated in TieBeNN.

| **Number assigned** | **Description** |
|:--------------------|:---------------:|
| 0 | IASP91 [(Kennett et al. 1995)](https://doi.org/10.1111/j.1365-246X.1995.tb03540.x) |
| 1 | AK135 [(Montagner & Kennett 1996)](https://doi.org/10.1111/j.1365-246X.1996.tb06548.x) |
| 2 | 2-layer-over-halfspace BGR velocity model [(Schlittenhardt 1999)](https://www.researchgate.net/profile/J-Schlittenhardt/publication/237600771_Regional_velocity_models_for_Germany_a_contribution_to_the_systematic_travel-time_calibration_of_the_international_monitoring_system/links/589dccbeaca272046aa92e2f/Regional-velocity-models-for-Germany-a-contribution-to-the-systematic-travel-time-calibration-of-the-international-monitoring-system.pdf) |
| 3 | WET, velocity model for Germany and adjacent regions |
| 4 | WB2012, Vogtland/ West Bohemia velocity model (1D, P+S), averaged from the 3D velocity models of [Ruzek & Horalek (2013)](https://doi.org/10.1093/gji/ggt295) |
| 5 | DEU, 3-layer model for Germany with a Moho depth of 28.5 km. It seems to be more appropriate for southern Germany. A vp/vs ratio of 1.68 was used to generate S-wave velocities |
| 6 | Crust1.0, for each latitude and longitude pair it produces up to 9 layers. Crust1.0 contains the topography at the epicenter. Python class created by Created by [J. Leeman and C. Ammon](https://github.com/jrleeman/Crust1.0) |
| 7 | Crust1.0 + AK135, the Crust1.0 model is used for the crustal structure and the lower continental structure from AK135 is merged immediately below |
| 8 | Crustal velocity model from [Küperkoch et al. (2018)](https://doi.org/10.1785/0120170365) for Insheim, (1D, P+S) |
| 9 | Regional model for LI (by J. Borns) |
|10 | Landau model (1D industry model, P+S) for the northern Upper Rhine Graben |
|11 | 1D model for Central Alps [(Diehl et al. 2021)](https://doi.org/10.1029/2021JB022155) |
|12 | 3D model for Central Alps (Diehl et al. 2021) _Note_: I am not sure I implemented this model correctly in NonLinLoc and 3D velocity models are currently deactivated |
|13 | 3D WEG model: shallow, high-resolution industry P-velocity model for the area around Rotenburg, Soehlingen, Soltau, Verden |
|14 | AlpsLocPS, 1D velocity models for the Greater Alpine Region [(Brazsus et al., 2024)](https://doi.org/10.1093/gji/ggae077) |
|15 | BENS, 1D P- and S-velocity models of the Northern Rhine Area [(Reamer & Hinzen 2004)](https://doi.org/10.1785/gssrl.75.6.713) |
|16 | ASZmod1, 1D P- and S-velocity models of the Albstadt Shear Zone [(Mader et al. 2021)](https://doi.org/10.5194/se-12-1389-2021) |
|17 | 3D P+S velocity models for the Upper Rhine Graben after [Lengline et al. (2023)](https://doi.org/10.1093/gji/ggad255) |
|18 | PO1, 1D P- and S-velocity model for the Novy-Kostel earthquake swarm region, constrained models [(Malek et al. 2023)](https://doi.org/10.1007/s00024-023-03250-w). Model defined down to 11 km depth |
|19 | PO1 (Malek et al. 2023) + WB2012 from 11 km depth |

When using `--vel_mode manual` the parameter `velmod` will be defined as the number assigned to the desired model. When the **Crust1.0** model is selected, the velocities in the layers at the epicenter are utilized to create a 1D velocity structure.

## Structure of NonLinLoc velocity structure files

Each layer uses one line in the text file defining the velocity models. These lines follow the NonLinLoc control file syntax:

````text
LAYER   DepthKm   VelP   GradVelP   VelS   GradVelS   Dens   GradDens
````

- `VelP`, `VelS` in km/s; `GradVelP`, `GradVelS` are the linear velocity gradients in km/s/km
- `Dens` in kg/m^3, `GradDens`: linear density gradient in kg/m^3/km

> :exclamation: **Note**
> The density is a carry over in the layer format from its original use in a waveform modelling code. It is not used in the NonLinLoc programs, so any convenient numerical value can be used.

## Add a new 1D velocity model

To add a new 1D velocity model in TieBeNN, you need to create a text file in the directory `<tiebenn_directory>/utils/velocity_models`. The file must named `v_<number_assigned>`
