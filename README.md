<p align='center'>
  <img src='figures/tiebenn_logo.png' />
</p>

## Description

:snail: **Here comes the documentation for TieBeNN!** :snail:

TieBeNN (**Tie**fen**Be**stimmung mittels **N**euronaler **N**etze) is an event-based wrapper which leverages a number of tools---some machine-learning based, some classic---to automatically produce phase picks for probabilistic hypocenter estimation of local events using [NonLinLoc](http://alomax.free.fr/nlloc/).

## Workflow

## Requirements

* **Python 3.9** or a later version.
* **SeisBench**, the popular seismology toolbox where the machine-learning models required by TieBeNN are stored.
* **NonLinLoc**, the set of programs written in C for probabilistic hypocenter estimation.
* **GMT**,
* **PyOcto**,
* **Pyrocko**,
* **GaMMA**,

## Installation

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

### Installing Python dependencies

Installing SeisBench will install most of Tiebenn's dependencies. You can install a pure-CPU version of SeisBench, in case it is necessary. For this, after activating the virtual environment by using the previously created alias, type:

```
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install seisbench
```

If you are working on a GPU machine, then you can skip the first line and directly install SeisBench with the second line. Then, you can proceed with the installation of the remaining dependencies:

```
pip install pygmt, pyocto, pyrocko, nllgrid
pip install git+https://github.com/wayneweiqiang/GaMMA.git 
```

Go to the desired directory where you wish to run Tiebenn and clone the repository:

```
git clone https://gitlab.szo.bgr.de/dzreorg/software/tiebenn.git
```

### Installing NonLinLoc and setting paths

### Install GMT

:snail: **...To be continued** :snail:

## Usage

### Input file

### Syntax

### Output

## To Do

This is a list of improvements I **would like** to implement in the code:

* Parallelize portions of the code
* Prepare the code documentation in Sphinx
* Fix NLL PDF visualization
* Add starplots to visualize the location quality based on diverse metrics

## Documentation

:pizza: :beer: Hungry for more detailed information? In-depth details about Tiebenn's functioning should have a future documentation (in the future).

## Authors and acknowledgment
C. Ramos

## Project status
Clearly under development.
