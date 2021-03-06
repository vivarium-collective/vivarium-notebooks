# Vivarium notebooks

## Installation
To run the code used in this work, you will need to install the following modules: `vivarium-core`, `vivarium-cobra`, 
`vivarium-bioscrape`, and `vivarium-multibody`.
These modules can be installed locally by executing the following command in the root directory:

> pip install -e ./

## Notebooks
Notebooks can be found under `notebooks/`. These include `Vivarium_interface_basics.ipynb` and `Multi-Paradigm-Composites.ipynb`.

## Python files
All Python files can be found under `bioscrape_cobra/`.
This includes Vivarium `Composers` for deterministic and stochastic versions of the Bioscrape/COBRA composite models 
called `bioscrape_cobra_deterministic.py` and `bioscrape_cobra_stochastic.py`. 
Simulation functions for running all of the examples in `Multi-Paradigm-Composites.ipynb` can be found in `simulate.py`.
This file also includes command-line run options for six different simulations with these names: 
[`deterministic`, `stochastic`, `deterministic_divide`, `stochastic_divide`, `deterministic_spatial`, `stochastic_spatial`].
These can be called simply with:

> python bioscrape_cobra/simulate --simulation_name

Parallelization can be triggered with the `-p` option: 

> python bioscrape_cobra/simulate --simulation_name -p

Saving the simulation output to a mongoDB database requires a running mongoDB instance, as described in the
[Vivarium docuentation](https://vivarium-core.readthedocs.io/en/latest/getting_started_dev.html). 
It can then be triggered with the `-d` option: 

> python bioscrape_cobra/simulate --simulation_name -d

