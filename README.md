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
* **GMT**
* **PyOcto**
* **Pyrocko**

Optional:

* **GaMMA**
* **nllgrid**


## Installation

I have tested TieBeNN in Linux Mint (and Lubuntu), so the instructions will use the Debian-based syntax.

It is highly recommended to use a *virtual environment* to install software with several requirements. We select the folder where we will install the virtual environment:

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
> Replace `<alias_name>` for something convenient. Then, reset the command line (or closing it and opening a new terminal also does it :laughing:). Next time you need to activate your virtual environment, just type `<alias_name>` in the command line.

### Installing dependencies

If you need a pure-CPU installation of SeisBench, you must do it manually. Installing SeisBench will install most of Tiebenn's dependencies. Thus, after activating the virtual environment by using the newly created alias, type:

```
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install seisbench
```

Go to the desired directory where you wish to run Tiebenn and clone the repository:

```
git clone https://gitlab.szo.bgr.de/dzreorg/software/tiebenn.git
```

:snail: **...To be continued** :snail:
srfrefer
## Usage

### Input file

### Syntax

### Output

## Documentation

:pizza::beer: Hungry for more detailed information? In-depth details about Tiebenn's functioning should be available (in the future) in the documentation.

## Authors and acknowledgment
C. Ramos

## Project status
Clearly under development.
