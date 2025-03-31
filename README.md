<p align='center'>
  <img src='figures/tiebenn_logo.png' />
</p>

## :memo: Description

TieBeNN (**Tie**fen**Be**stimmung mittels **N**euronaler **N**etze) is an event-based wrapper which leverages a number of tools---some machine-learning based, some classic---to automatically produce phase picks for probabilistic hypocenter estimation of local events using [NonLinLoc](http://alomax.free.fr/nlloc/).

## :gear: Workflow

Using as input the coordinates (latitude and longitude) and the time of a local event, the seismic event location with TieBeNN goes through the following stages:

1. **Waveform data fetching**: It produces a catalog of stations around the epicenter. Then, waveform data in miniSEED format are retrieved using an ObsPy client on FDSN servers, or from a SDS directory structure.
1. **Waveform data preprocessing**: empty channels are removed, masked channels are split, data are detrended and bandpass filtered. Optionally, stations within 100 km from the epicenter are denoised using the model DeepDenoiser [(Zhu et al. 2019)](https://arxiv.org/abs/1811.02695).
1. **Phase picking**: P- and S-phases (first arrivals) go through the phase-picking models (either EQTransformer and PhaseNet, as available in the [SeisBench](https://github.com/seisbench/seisbench) toolbox) for phase detection.
1. **Phase association**: Phase picks go through phase association to discard possible false detections. The associators available in TieBeNN are [PyOcto](https://github.com/yetinam/pyocto) and [GaMMA](https://github.com/AI4EPS/GaMMA).
1. **Outputs exported**: Event-associated phase picks on each station are exported in a CSV file, including station coordinates, phase arrival time, probability and signal-to-noise ratio. Input files requested by NonLinLoc are generated. Optionally, figures are generated of waveforms and the detected picks at each station, as well as phase association figures.
1. **Probabilistic hypocenter estimation**: The generated files are used for hypocenter estimation with NonLinLoc. Optionally, figures are generated: :one: epicenter and stations with picks in map :two: waveforms with picks, sorted by epicentral distances :three: confidence ellipsoid of event location.
1. **Location quality assessment**: Location metrics are gathered/calculated to generate the Location Quality Score metric, as well as a figure to visualize it. :memo: **A description of this metric will be available in a manuscript, currently in preparation** :memo:

> :point_right: **Note**: When TieBeNN starts, it will go on in a loop until the minimum requested detections within a given epicentral distance are obtained. If this is not the case (there could be several reasons for this), the epicentral distance for station waveform retrieval is gradually increased and the process is repeated. If not enough phase picks are detected within 200 km, the loop breaks and TieBeNN end the run reporting an unsuccessful event location.

## :white_check_mark: Requirements

* **Python 3.9** or a later version.
* **SeisBench**, the popular seismology toolbox where the machine-learning models required by TieBeNN are stored.
* **NonLinLoc**, the set of programs written in C for probabilistic hypocenter estimation.
* **GMT**,
* **PyOcto**,
* **Pyrocko**,
* **GaMMA**,

## :hammer_and_wrench: Installation

I have tested TieBeNN in Linux Mint (and Lubuntu), so the instructions will use the Debian-based syntax.

It is highly recommended to use a *virtual environment* to install software with several requirements. We select the folder where we will install the virtual environment (do not forget to replace the names between angle brackets for what you choose):

```
python3 -m venv <path_to_virtual_environment>/<venv_tiebenn>
```

Then, we must activate the virtual environment:

```
source <path_to_virtual_environment>/<venv_tiebenn>/bin/activate
```
> :bulb: **TIP**
>
> You can create a shortcut in your _.bashrc_ file to quickly access the virtual environment in future runs. Open _.bashrc_ with your favorite editor (e.g. `nano ~/.bashrc`) and add the following line at the end of the file:
>
> `alias <alias_name>='source <path_to_virtual_environment>/<venv_tiebenn>/bin/activate'`
>
> Replace `<alias_name>` for something convenient. Then, reset the command line:
> `exec bash`

> for changes to take effect. To activate your virtual environment, just type `<alias_name>` in the command line.

### :hammer_and_wrench: Installing Python dependencies

Installing SeisBench will install most of Tiebenn's dependencies. You can install a pure-CPU version of SeisBench, in case it is necessary. For this, after activating the virtual environment by using the previously created alias, type:

```
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install seisbench
```

If you are working on a GPU machine, then you can skip the first line and directly install SeisBench with the second line. Then, you can proceed with the installation of the remaining dependencies:

```
pip install pygmt pyocto pyrocko nllgrid
pip install git+https://github.com/wayneweiqiang/GaMMA.git 
```

Go to the desired directory where you wish to run Tiebenn and clone the repository:

```
git clone https://gitlab.szo.bgr.de/dzreorg/software/tiebenn.git
```

### :hammer_and_wrench: Installing NonLinLoc and setting paths



### :hammer_and_wrench: Install GMT

:snail: **...To be continued** :snail:

Pygmt dice que necesita >=6.4.0, yo tengo 6.6.0, gpu01f tiene 6.0.0

## :test_tube: Usage

### :inbox_tray: Input file

### :receipt: Syntax

### :outbox_tray: Output

## :wrench: To Do

This is a list of improvements I **would like** to implement in the code:

* Prepare the code documentation in Sphinx
* Prepare a containerized version of the code

## :books: Documentation

:pizza: :beer: Hungry for more detailed information? In-depth details about Tiebenn's functioning should have a future documentation (in the future).

## :book: References

A software report is currently in preparation.

## :brain: Authors and acknowledgment
C. Ramos

## :construction: Project status :construction:
Clearly under development.

## :link: Useful links

[DeepDenoiser example] (https://colab.research.google.com/github/seisbench/seisbench/blob/main/examples/02b_deep_denoiser.ipynb)
[NonLinLoc GitHub] https://github.com/ut-beg-texnet/NonLinLoc
