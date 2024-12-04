<p align='center'>
  <img src='figures/tiebenn_logo.png' />
</p>

## Description

:snail: **Here comes the documentation for TieBeNN!** :snail:

TieBeNN (**Tie**fen**Be**stimmung mittels **N**euronaler **N**etzwerke) is an event-based wrapper which leverages a number of tools---some machine-learning-based, some classic---to automatically produce phase picks for probabilistic hypocenter estimation of local events using [NonLinLoc](http://alomax.free.fr/nlloc/).

## Requirements

* **Python 3.9** or a later version.
* **SeisBench**, the popular seismology toolbox where the machine-learning models required by TieBeNN are stored.

## Installation

I have tested TieBeNN in Linux Mint (and Lubuntu), so the instructions will use the Debian-based syntax.

It is highly recommended to use a *virtual environment* to install software with several requirements. We select the folder where we will install the virtual environment:

```
python3 -m venv <path_to_virtual_environment>/venv_tiebenn
```

Then, we must activate the virtual environment:

```
source <path_to_virtual_environment>/venv_tiebenn/bin/activate
```
> :bulb: **TIP**
>
> You can create a shortcut in your _.bashrc_ file to quickly access the virtual environment in future runs. Open _.bashrc_ with your favorite editor (e.g. `nano ~/.bashrc`) and add the following line at the end of the file:
>
> `alias <alias_name>='source <path_to_virtual_environment>/venv_tiebenn/bin/activate'`
>
> Replace `<alias_name>` for something convenient. Then, reset the command line (or closing it also does it :laughing:). Next time you need to activate your virtual environment, just type `<alias_name>` in the command line.

:snail: **...To be continued** :snail:

## Authors and acknowledgment
C. Ramos

## Project status
Clearly under development.
