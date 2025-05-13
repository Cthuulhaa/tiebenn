<p align="center">
  <img src="figures/tiebenn_logo.svg" width="600"/>
</p>

## :memo: Description

TieBeNN (**Tie**fen**Be**stimmung mittels **N**euronaler **N**etze) is an event-based wrapper that leverages several tools—some machine-learning-based, some traditional—to automatically generate phase picks for probabilistic hypocenter estimation of local events using [NonLinLoc](http://alomax.free.fr/nlloc/).

## :gear: Workflow

Using the coordinates (latitude and longitude) and the UTC time of a local event as input, TieBeNN processes seismic event locations through the following stages:

1. **Waveform data fetching**: A catalog of stations around the epicenter is produced. Then, waveform data in miniSEED format are retrieved using an ObsPy client on FDSN servers or from a SDS directory structure.

1. **Waveform data preprocessing**: Empty channels are removed, masked channels are split, and data are detrended and bandpass filtered. Optionally, stations within 100 km from the epicenter are denoised using the DeepDenoiser model [(Zhu et al. 2019)](https://arxiv.org/abs/1811.02695).

1. **Phase picking**: P- and S-phases (first arrivals) are detected using phase-picking models (either EQTransformer or PhaseNet, available in the [SeisBench](https://github.com/seisbench/seisbench) toolbox).

1. **Phase association**: Detected picks are passed through phase associators to discard possible false detections. TieBeNN supports [PyOcto](https://github.com/yetinam/pyocto) and [GaMMA](https://github.com/AI4EPS/GaMMA).

1. **Export outputs**: Event-associated phase picks are exported to CSV files, including station coordinates, arrival times, pick probabilities, and signal-to-noise ratios. NonLinLoc input files are generated. Optionally, waveform and pick figures per station, as well as phase association plots, are produced.

1. **Probabilistic hypocenter estimation**: The generated files are used by NonLinLoc for hypocenter estimation. Optionally, figures are created: :one: epicenter and stations on a map; :two: waveforms with picks, sorted by epicentral distance; :three: location confidence ellipsoid.

1. **Location quality assessment**: Location metrics are gathered to compute the Location Quality Score (LQS) and to generate a visualization. :memo: **A description of this metric should be available in a manuscript, currently in preparation** :memo:

> :point_right: **Note**:
>
> TieBeNN loops through this process until the minimum required detections within a given epicentral distance are obtained. If not, the search radius is gradually expanded. If phase picks are insufficient within 200 km, the run ends with an unsuccessful event location.

## :white_check_mark: Requirements

* **Python 3.9** or later (successfully tested with Python 3.10 and 3.12).
* **SeisBench**, the ML model toolbox used for phase picking and denoising.
* **NonLinLoc**, a suite of C programs for probabilistic hypocenter estimation.
* **PyOcto**, phase associator after [Münchmeyer (2024)](https://seismica.library.mcgill.ca/article/view/1130).
* **Pyrocko**, open-source seismology toolbox and library.
* **GaMMA**, phase associator after [Zhu et al. (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023249).
* **Cartopy**, for map generation.

### Optional

* **NLLGrid**, a Python class for handling NonLinLoc grid files. Hosted [here](https://github.com/claudiodsf/nllgrid). Useful if you would eventually like to try 3D NonLinLoc grids in Python.
* **Sphinx**, in case you want to generate TieBeNN's documentation.

## :hammer_and_wrench: Installation

Tested on Linux Mint and Lubuntu (Debian-based syntax used here).

### :hammer_and_wrench: Create a virtual environment

It is highly recommended to use a *virtual environment* to install software with several requirements. We move to the folder where we will install the virtual environment (do not forget to replace the names between angle brackets). The created environment can be then activated using `source`:

```bash
python3 -m venv <path_to_virtual_environment>/<venv_tiebenn>
source <path_to_virtual_environment>/<venv_tiebenn>/bin/activate
```

> :bulb: **TIP**
>
> You can add an alias by including the following line at the end of your `~/.bashrc` for quick access:
>
> ```bash
> alias <alias_name>='source <path_to_virtual_environment>/<venv_tiebenn>/bin/activate'
> ```
> Save changes. Then reload:
> ```bash
> exec bash
> ```

### :hammer_and_wrench: Installing TieBeNN and its Python dependencies

In case you plan to work on a CPU machine, after activating the virtual environment by using the previously created alias, type:

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

If you are working on a GPU machine, then you can skip this line. Then, you can proceed to clone the TieBeNN repository and install it:

```bash
git clone https://github.com/Cthuulhaa/tiebenn
cd tiebenn

pip install .
pip install -r optional.txt
```

> :point_right: **Note**
>
> Your are free to modify or comment-out those packages you will not need.

### :hammer_and_wrench: Installing NonLinLoc and setting paths

_(Follow these steps in case NonLinLoc is not installed)_ NonLinLoc must be individually compiled to make sure it is compatible with the machine where TieBeNN will be running. First we will create a convenient directory to install NonLinLoc and within it, we clone the NonLinLoc repository:

```bash
mkdir -p <some_convenient_directory>
cd <some_convenient_directory>
git clone https://github.com/ut-beg-texnet/NonLinLoc.git
cd NonLinLoc/src
mkdir bin
cmake .
make
cd bin
cp Vel2Grid* Grid2* NLLoc ../../../
```

> :exclamation: **Important**
>
> Do not use NonLinLoc's latest release directly, as it might contain unresolved bugs, whose fix are still unreleased.

Set NonLinLoc in your `PATH`:

```bash
echo 'export PATH=${PATH}:<some_convenient_directory>' >> ~/.bashrc
exec bash
```

## :test_tube: Usage

This section shows an example which should make TieBeNN's usage clear.

### :inbox_tray: Input file

Input file syntax:
```text
YYYY-MM-DDTHH:MM:SS  latitude  longitude
```
Example:
```text
2024-12-24T01:55:04 51.337  12.548
```

### :receipt: Syntax

With the virtual environment activated (see above), TieBeNN follows this syntax:

```bash
tiebenn --event_file <EventFile> --max_epic_dist <MaxEpDist> --picker <Picker> --client <Client> --sds_dir <SDSDir> --min_detections <MinDetections> --plots <Plots> --vel_mode <VelMode> --velmod <VelMod> --ph_assoc <PhaseAssoc> --denoise <Denoise> --mult_windows <MultiWindows> --secs_before <SecsBef>
```

| Parameter | Description |
|:----------|:-----------:|
| **`--event_file`** | Full path to input file with preliminary epicenter (latitude, longitude) and UTC datetime |
| **`--max_epic_dist`** | Maximum epicentral distance (in km) for stations on which phase picks will be detected |
| **`--picker`** | Select model for phase picking. **EQTransformer** can be defined as `sb_eqt`, `seisbench_eqt`, `seisbench_eqtransformer`, `sb_eqtransformer` and **PhaseNet** as `sb_pn`, `seisbench_pn`, `sb_phasenet`, `seisbench_phasenet`. Not case sensitive |
| **`--client`** | If set to `SDS`, it will access a directory with SeisComp3 structure. It will try to fetch the available waveforms from the stations in the station list within `MaxEpDist` km. Afterwards, it will try to fetch stations using FDSN clients to access their services. If set to `FDSN`, it skips the search for a SDS directory |
| **`--sds_dir`** | A string with the full path to the SeisComp3 directory. This parameter must be defined if `Client` is set to `SDS` |
| **`--min_detections`** | Minimum amount of stations on which P- or S- phase picks must be detected for the detection loop to end |
| **`--plots`** | If set to True, it will plot the waveforms recorded on each stations with at least one phase detection. It will also plot all the phase picks associated to the event sorted by epicentral distance, as well as plots of the locations: epicenter and stations with detections on a map, waveforms with phase picks sorted by epicentral distance, and confidence ellipsoid of event location |
| **`--vel_mode`** | This parameter decides how to choose a seismic velocity model for the event location. Options are `automatic` (`automatic`, `auto` or `a`) for choosing a local velocity model based on the epicenter (if available, otherwise the layer-over-halfspace model used at EdB is used as default) and `manual` (`manual`, `man`, `m`) for choosing the velocity model manually. In this case, the parameter `VelMod` must be defined |
| **`--velmod`** | The number corresponding to the model which will be used for hypocenter location with NonLinLoc. [See the full list here](velocity_models.md) _Note_: 3D velocity models have been tested, although I have still not found a region in Germany where using a 3D velocity model instead of a dedicated, local 1D velocity model results in a dramatic improvement in the event location quality and is worth the extra travel-time calculation time. The implementation of 3D velocity models for real-time event location is, at least momentarily, beyond the scope of this repository |
| **`--ph_assoc`** | Phase associator. Options are `PyOcto` (`pyocto` or `p`) and `GaMMA` (`gamma` or `g`). Not case sensitive |
| **`--denoise`** | Boolean parameter. If true, the DeepDenoiser model will be applied on the waveforms of stations within 100 km in epicentral distance |
| **`--mult_windows`** | Boolean parameter. If true, it makes the phase picker to look for P- or S-waves in moving windows, which helps to address the prediction inconsistency inherent to machine-learning-based models |
| **`--secs_before`** | If **`--mult_windows`** is _False_, then you can use this parameter to set the starttime of the retrieved waveforms (in seconds before the event time). Default is 0 seconds |

> :point_right: **Note**:
>
> Commands related to the use of 3D NonLinLoc velocity models is deactivated, but they are still present in the code in case someone wants to pick its development up.

To locate the event with example UTC datetime and coordinates specified above, we type in the terminal (remember to activate the virtual environment created above):

```bash
tiebenn --event_file <full_path_to_example_event> --max_epic_dist 150 --picker SeisBench_PhaseNet --client FDSN --min_detections 3 --plots True --vel_mode auto --ph_assoc PyOcto --denoise True --mult_windows True
```

### :outbox_tray: Output

If the location was successful, you’ll see something like:
```
##############################################################
NonLinLoc: Location completed.
Origin time: 24-12-2024_01:55:4.517789
             Lat:51.328045 Long:12.543796 Depth: 22.685547 km
--------------------------------------------------------------
Location quality:
RMS:0.285177 Number of phases: 130 Gap: 35.032 Distance from hypocenter to nearest station: 12.077 km
Ellipsoid semi-major axis: 1.030196e+00
##############################################################
```

With `--plots True`, the following figures are also generated:

* Phase picks on waveforms:
  <p align="center">
    <img src="figures/example_picks.svg" width="500"/>
  </p>

* Associated phases:
  <p align="center">
    <img src="figures/example_phassoc.svg" width="500"/>
  </p>

* Map: epicenter and stations:
  <p align="center">
    <img src="figures/example_epicenter_map.png" width="500"/>
  </p>

* Location confidence ellipsoid:
  <p align="center">
    <img src="figures/example_ellipsoid.png" width="500"/>
  </p>

* Location Quality Score:
  <p align="center">
    <img src="figures/example_LQS.svg" width="400"/>
  </p>

> :point_right: **Ideas welcome!** Please submit feature suggestions via new issues.

## :books: Documentation

:pizza: :beer: Hungry for more detailed information? The full documentation is [here](https://tiebenn.readthedocs.io/en/latest/)

## :book: References

If you use TieBeNN in your research, please cite:

> TieBeNN: A neural network-based tool for automatic focal depth estimation
>
> Ramos, C. (2025) _TieBeNN v0.2.0_ [Software]. Zenodo. [https://doi.org/10.5281/zenodo.15384316](https://doi.org/10.5281/zenodo.15384316)

A detailed paper/software report on TieBeNN is currently in preparation.

## :brain: Authors and Acknowledgment

C. Ramos (maintainer)

## :construction: Project Status

Clearly under development.

## :link: Useful Links

- [DeepDenoiser example](https://colab.research.google.com/github/seisbench/seisbench/blob/main/examples/02b_deep_denoiser.ipynb)
- [NonLinLoc GitHub](https://github.com/ut-beg-texnet/NonLinLoc)
- [Pyrocko applications](https://pyrocko.org/)

## :balance_scale: License

This project is released under the **GNU General Public License v3.0 (GPLv3)**.
In short, it's open-source and free to use, modify, and redistribute — but any modified version that is shared must also remain open under the same license.

See the full [LICENSE](LICENSE) file for legal details.

### :mag: TL;DR (not legally binding!)

- :white_check_mark: You can use, modify, and share this code freely.
- :scroll: If you share a modified version, you **must** also share the source code.
- :brain: The license ensures that TieBeNN and its derivatives stay open-source.
- :x: No warranties — the software is provided *as is*.
