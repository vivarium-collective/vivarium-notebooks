# Vivarium notebooks

This repository contains the simulations described in: 
[Agmon, E., Spangler, R.K., Skalnik, C.J., Poole, W., Morrison, J.H., Peirce, S.M., and Covert, M.W. (2021). 
Vivarium: an interface and engine for integrative multi-scale modeling in computational biology. *Bioinformatics*.](
https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btac049/6522109)

[Vivarium Core](https://github.com/vivarium-collective/vivarium-core) is required to start your own Vivarium project. 
See the [documentation](https://vivarium-core.readthedocs.io/en/latest/getting_started.html) for more details,
and the [Vivarium Collective](https://vivarium-collective.github.io) for available modules.


## Notebooks
Notebooks can be found under `notebooks/`. These include:

 * [Lac Operon CRN](https://nbviewer.org/github/vivarium-collective/vivarium-notebooks/blob/main/notebooks/Lac_Operon_CRN.ipynb)
 * [Multi-Paradigm Composites](https://nbviewer.org/github/vivarium-collective/vivarium-notebooks/blob/main/notebooks/Multi-Paradigm-Composites.ipynb)
 * [Vivarium Interface Basics](https://nbviewer.org/github/vivarium-collective/vivarium-core/blob/master/notebooks/Vivarium_interface_basics.ipynb)

## Setup
To run the code used in this work, you will need to install the following modules: `vivarium-core`, `vivarium-cobra`, 
`vivarium-bioscrape`, and `vivarium-multibody`.
These modules can be installed locally by executing the following command in the root directory:

1. Clone the repository.

   ```console
   $ git clone https://github.com/vivarium-collective/vivarium-notebooks.git
   ```

   Alternatively, you can simply extract an archive of the source code
   if you got the code in that format.

2. (recommended but optional) Setup a Python virtual environment.
3. Install Python dependencies

   ```console
   $ pip install -r requirements.txt
   ```

   Note that you have to install numpy first because the `setup.py`
   script of one of our dependencies requires it.

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
[Vivarium documentation](https://vivarium-core.readthedocs.io/en/latest/getting_started_dev.html). 
It can then be triggered with the `-d` option: 

> python bioscrape_cobra/simulate --simulation_name -d

## Tests
Tests are performed with pytest. Simply call the following to ensure everything is working properly:
```
$ pytest
```

To run only the fast tests:
```
$ pytest -m 'not slow'
```

