"""
==============================================
Simulation helper functions for BioscrapeCOBRA
==============================================

Includes functions for configuring, running, and plotting all experiments reported in the paper:
    "Vivarium: an interface and engine for integrative multiscale modeling in computational biology"

These are used in the supplementary Jupyter notebook "Multi-Paradigm-Composites.ipynb".

This file also includes pre-configured experiments that can be triggered from the command line:
    * 1 : deterministic single cell
    * 2 : stochastic single cell
    * 3 : deterministic with cell division
    * 4 : stochastic with cell division
    * 5 : deterministic cells in spatial environment
    * 6 : stochastic cells in spatial environment

These can be triggered from the command line by entering the simulation number:
```
$ python bioscrape_cobra/simulate.py -sim_number
```

parallel and mongoDB database for saved output can be called with:
```
$ python bioscrape_cobra/simulate.py [-database, -d] [-parallel, -p]
```

see all options with:
```
python bioscrape_cobra/simulate.py -h
```

"""
import os
import argparse
import copy
import random
import time as clock
import pytest
from tqdm import tqdm

# vivarium imports
from vivarium.core.engine import Engine, timestamp
from vivarium.library.units import units
from vivarium.library.dict_utils import deep_merge
from vivarium.core.composition import simulate_composer, simulate_composite
from vivarium.core.composer import Composite
from vivarium.library.units import remove_units
from vivarium.core.serialize import deserialize_value

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import (
    Lattice, make_lattice_config)
from vivarium_multibody.composites.grow_divide import GrowDivide
from vivarium_multibody.processes.diffusion_field import DiffusionField

# cobra imports
from vivarium_cobra.processes.cobra_fba import (
    COBRA_FBA, 
    get_iAF1260b_config, 
)
from vivarium_cobra.composites.cobra_composite import CobraComposite
from vivarium_cobra.library.lattice_utils import get_bin_volume

# Bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# local import
from bioscrape_cobra.bioscrape_cobra_stochastic import (
    BioscrapeCOBRAstochastic, SBML_FILE_STOCHASTIC, COBRA_TIMESTEP, BIOSCRAPE_TIMESTEP)
from bioscrape_cobra.bioscrape_cobra_stochastic import GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL
from bioscrape_cobra.bioscrape_cobra_deterministic import BioscrapeCOBRAdeterministic, SBML_FILE_DETERMINISTIC

# plotting
from vivarium.plots.topology import plot_topology
from vivarium.plots.simulation_output import _save_fig_to_dir as save_fig_to_dir
from bioscrape_cobra.plot import (
    plot_multigen, plot_single, plot_fields_tags, plot_fields_snapshots, config_embedded_bioscrape_cobra_topology)

# get the sbml files at the correct path
dirname = os.path.dirname(__file__)
sbml_deterministic_file = os.path.join(dirname, SBML_FILE_DETERMINISTIC)
sbml_stochastic_file = os.path.join(dirname, SBML_FILE_STOCHASTIC)


# default environment variables, which can be varied by the function `simulate_bioscrape_cobra`
DEFAULT_DIVIDE_THRESHOLD = 2000 * units.fg
INITIAL_GLC = 10  # mmolar
INITIAL_LAC = 20  # mmolar
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 2

# plotting definitions
plot_variables_list = [
    ('species', 'rna_M'),
    ('species', 'monomer_betaGal'),
    ('species', 'protein_betaGal'),
    ('species', 'protein_Lactose_Permease'),
    ('species', GLUCOSE_EXTERNAL),
    ('species', LACTOSE_EXTERNAL),
    ('flux_bounds', 'EX_glc__D_e'),
    ('flux_bounds', 'EX_lac__D_e'),
    ('boundary', 'no_units', 'biomass'),
]
plot_variables_list_deterministic = [
    ('rates', 'k_dilution__',)]
plot_variables_list_deterministic.extend(plot_variables_list)
plot_variables_list_stochastic = []
plot_variables_list_stochastic.extend(plot_variables_list)


# Simulates Cobra Model on its Own
def simulate_cobra(
    total_time=100,
    initial_state=None,
):
    # get the configuration for the iAF1260b BiGG model
    iAF1260b_config = get_iAF1260b_config()
    iAF1260b_config.update({'time_step': COBRA_TIMESTEP})

    iAF1260b_config['name'] = 'cobra' #Rename the process
    dynamic_fba = COBRA_FBA(iAF1260b_config)

    cobra_composite = Composite({
        'processes': dynamic_fba, 
        'topology': dynamic_fba.generate_topology()
        })

    # get the initial state
    if initial_state is None:
        initial_state = dynamic_fba.initial_state({})

    # run simulation
    cobra_sim_settings = {
        'initial_state': initial_state,
        'total_time': total_time}
    cobra_timeseries = simulate_composer(dynamic_fba, cobra_sim_settings)

    return cobra_timeseries, dynamic_fba


# Simulates Cobra Model in a Composite with some Derivers
def simulate_cobra_composite(
    total_time=100,
    initial_state=None,
):

    # get the configuration for the iAF1260b BiGG model
    iAF1260b_config = get_iAF1260b_config()
    iAF1260b_config.update({'time_step': COBRA_TIMESTEP})

    #Place the Process into a Composite level dictionary
    composite_config = {'cobra': iAF1260b_config}
    cobra_composite = CobraComposite(composite_config)

    # get the initial state
    if initial_state is None:
        initial_state = cobra_composite.initial_state({})

    # run simulation
    cobra_sim_settings = {
        'initial_state': initial_state,
        'total_time': total_time}
    cobra_timeseries = simulate_composer(cobra_composite, cobra_sim_settings)

    return cobra_timeseries, cobra_composite


def simulate_bioscrape(
        stochastic=False,
        initial_glucose=None,
        initial_lactose=None,
        initial_state=None,
        total_time=100,
        initial_volume=1.0,
        sbml_file=None,
):
    
    #create configs
    if not stochastic:
        bioscrape_config = {
            'sbml_file': sbml_file or sbml_deterministic_file,
            'stochastic': False,
            'initial_volume': initial_volume,
            'internal_dt': BIOSCRAPE_TIMESTEP/100,
            'time_step': BIOSCRAPE_TIMESTEP}
    else:
        bioscrape_config = {
            'sbml_file': sbml_file or sbml_stochastic_file,
            'stochastic': True,
            'safe_mode': False,
            'initial_volume': initial_volume,
            'internal_dt': BIOSCRAPE_TIMESTEP/100,
            'time_step': BIOSCRAPE_TIMESTEP
        }
    
    bioscrape_process = Bioscrape(bioscrape_config)


    #Create an Empty Composite for Plotting Purposes
    bioscrape_composite = Composite({
        'processes': bioscrape_process.generate_processes(), 
        'topology': bioscrape_process.generate_topology()
        })

    # get the initial state
    if initial_state is None:
        initial_state = bioscrape_process.initial_state({})

    if initial_glucose is not None:
        initial_state['species']['Glucose_external'] = initial_glucose

    if initial_lactose is not None:
        initial_state['species']['Lactose_external'] = initial_lactose

    # run simulation
    bioscrape_sim_settings = {
        'initial_state': initial_state,
        'total_time': total_time
        }

    bioscrape_timeseries = simulate_composer(bioscrape_process, bioscrape_sim_settings)

    return bioscrape_timeseries, bioscrape_composite


# Simulate a system of cells that grow and divide in a well mixed spatial environment
def simulate_grow_divide(
        total_time=100,
        growth_rate=0.03,
        initial_state=None,
        growth_noise=10**-6
):
    # config
    growth_config = {
        'default_growth_rate': growth_rate,
        'default_growth_noise': growth_noise,
        'agents_path': ('agents',)}
    grow_divide_composer = GrowDivide({
        'agent_id': '0', 'growth' : growth_config})

    if initial_state is None:
        initial_state = {
            'agents': {
                '0': {
                    'global': {
                        'mass': 1000 * units.femtogram}
                }}}

    grow_divide_sim_settings = {
        'total_time': total_time,
        'initial_state': initial_state,
        'outer_path': ('agents', '0'),
        'return_raw_data': True}
    grow_divide_composite = grow_divide_composer.generate(path=('agents', '0'))

    grow_divide_data = simulate_composite(grow_divide_composite, grow_divide_sim_settings)
    grow_divide_data = deserialize_value(grow_divide_data)
    grow_divide_data = remove_units(grow_divide_data)

    #returns the data, the initial composite, and the final composite
    return grow_divide_data, grow_divide_composer.generate(path=('agents', '0')), grow_divide_composite


def simulate_diffusion(
        total_time=100,
        diffusion_rate=0.001,
        initial_state={},
        bins=[10, 10],
        bounds=[10, 10]
):
    # configure
    

    config = {
        'n_bins':bins,
        'bounds':bounds,
        'diffusion':diffusion_rate,
        'initial_state':{'glc':initial_state}
    }

    diffusion_process = DiffusionField(config)

    diffusion_composer = Composite({
        'processes': diffusion_process.generate_processes(), 
        'topology': diffusion_process.generate_topology()
        })

    diffusion_state = diffusion_process.initial_state({})

    diffusion_sim_settings = {
        'total_time': total_time,
        'return_raw_data': True,
        'initial_state': diffusion_state,
    }

    diffusion_data = simulate_composer(diffusion_process, diffusion_sim_settings)

    #Add empty agents to reuse plotting functionality
    for t in diffusion_data:
        diffusion_data[t]['agents'] = {}

    return diffusion_data, diffusion_composer


def get_lattice_grow_divide_composite(
        diffusion_rate=0.001,
        initial_concentration={},
        bins=[10, 10],
        bounds=[10, 10],
        growth_rate=0.03,
        growth_noise=10**-6,
        depth=10
):
    lattice_config = make_lattice_config(
            bounds=bounds,
            n_bins=bins,
            concentrations={'glc':initial_concentration},
            diffusion=diffusion_rate,
            depth=depth)

    lattice_composer = Lattice(lattice_config)
    lattice_composite = lattice_composer.generate()


    growth_config = {'default_growth_rate': growth_rate, 'default_growth_noise': growth_noise}
    grow_divide_composer = GrowDivide({'agent_id': '0', 'growth': growth_config})

    agent_id = '0'
    lattice_grow_divide_composite = grow_divide_composer.generate(path=('agents', agent_id))
    lattice_grow_divide_composite.merge(composite=lattice_composite)

    return lattice_grow_divide_composite


def simulate_grow_divide_lattice(
        lattice_grow_divide_composite,
        total_time=100,
        initial_state=None,
):

    agent_id = '0'
    if initial_state is None:
        initial_state = {
            'agents': {
                agent_id: {
                    'global': {
                        'mass': 1000 * units.femtogram}
                }}}

    sim_settings = {
        'total_time': total_time,
        'return_raw_data': True,
        'initial_state': initial_state,
        'return_raw_data': True}
    lattice_grow_divide_data = simulate_composite(
        lattice_grow_divide_composite, sim_settings)

    return lattice_grow_divide_data


# helper functions
def get_bioscrape_cobra_config(
        spatial=False,
        division=False,
        divide_threshold=DEFAULT_DIVIDE_THRESHOLD,
        external_volume=None,
        sbml_file=None,
        parallel=False,
):
    """ create a generic config dict for bioscrape_cobra composers """
    config = {
        '_parallel': parallel,
        **({'local_fields': {'bin_volume': external_volume}}
           if external_volume is not None else {}),
        **({'sbml_file': sbml_file} if sbml_file is not None else {}),
    }

    if spatial or division:
        config.update({
            'divide_on': True,
            'agents_path': ('..', '..', 'agents',),
            'fields_path': ('..', '..', 'fields',),
            'dimensions_path': ('..', '..', 'dimensions',),
            'divide_condition': {
                'threshold': divide_threshold},
        })

    if spatial:
        config.update({
            'fields_on': True})

    return config


def get_initial_state(initial_states, n=None):
    """ return an initial state from the initial_states argument """
    state_copy = copy.deepcopy(initial_states)
    if isinstance(state_copy, list):
        if n is not None and len(state_copy) > n:
            return state_copy[n]
        else:
            return random.choice(state_copy)
    elif isinstance(state_copy, dict):
        return state_copy
    else:
        return {}


def simulate_bioscrape_cobra(
        division=False,
        stochastic=False,
        spatial=False,
        initial_glucose=1e1,
        initial_lactose=1e1,
        initial_agent_states=None,
        bounds=None,
        n_bins=None,
        depth=DEPTH,
        diffusion_rate=1e-1,
        jitter_force=1e-4,
        divide_threshold=2000 * units.fg,
        external_volume=None,
        n_agents=1,
        halt_threshold=100,
        total_time=100,
        sbml_file=None,
        emitter='timeseries',
        output_type=None,
        parallel=False,
):
    """ Main simulation function for BioscrapeCOBRA

    Args:
        * division:
        * stochastic:
        * spatial:
        * initial_glucose:
        * initial_lactose:
        * initial_agent_states:
        * bounds:
        * n_bins:
        * depth:
        * diffusion_rate:
        * divide_threshold:
        * external_volume:
        * n_agents:
        * halt_threshold:
        * total_time:
        * sbml_file:
        * emitter:
        * output_type:
        * parallel:
    """

    if n_bins is None:
        n_bins = NBINS
    if bounds is None:
        bounds = BOUNDS

    # get the bin volume based upon the lattice
    bin_volume = (external_volume or get_bin_volume(n_bins, bounds, depth)) * units.L
    agent_state = {
        # set initial flux values based on COBRA defaults.
        'flux_bounds': {
            'EX_lac__D_e': 0.0,
            'EX_glc__D_e': 0.099195
        },
        # field_counts_deriver needs the bin volume
        'boundary': {
            'bin_volume': bin_volume,
            **({'external': {
                GLUCOSE_EXTERNAL: initial_glucose,
                LACTOSE_EXTERNAL: initial_lactose}
               } if spatial else {}),
        },
    }

    # make the BioscrapeCOBRA config
    biocobra_config = get_bioscrape_cobra_config(
        spatial=spatial,
        division=division,
        divide_threshold=divide_threshold,
        external_volume=bin_volume,
        sbml_file=sbml_file,
        parallel=parallel)

    # get the BioscrapeCOBRA composer -- either stochastic or deterministic
    if stochastic:
        biocobra_composer = BioscrapeCOBRAstochastic(biocobra_config)
    else:
        biocobra_composer = BioscrapeCOBRAdeterministic(biocobra_config)

    # initialize the composite
    biocobra_composite = Composite({
        'processes': {},
        'topology': {}})

    if spatial:
        # make n_agents dividing agents and combine them with a Lattice environment

        agents_initial = {'agents': {}}
        for n in range(n_agents):
            agent_id = str(n)
            config = {'agent_id': agent_id}
            agent = biocobra_composer.generate(config=config)
            biocobra_composite.merge(
                composite=agent,
                path=('agents', agent_id))

            # initial state for agent
            agents_initial['agents'][agent_id] = get_initial_state(
                initial_agent_states, n=n)
            agents_initial['agents'][agent_id] = deep_merge(
                agents_initial['agents'][agent_id], agent_state)

        # create a second initial composite for plotting
        initial_composite = biocobra_composer.generate(
            path=('agents', '0'),
            config={'agent_id': '0'})

        # make a lattice composite for the environment
        field_concentrations = {
            GLUCOSE_EXTERNAL: initial_glucose,
            LACTOSE_EXTERNAL: initial_lactose}
        lattice_config = make_lattice_config(
            bounds=bounds,
            n_bins=n_bins,
            depth=depth,
            concentrations=field_concentrations,
            diffusion=diffusion_rate,
            jitter_force=jitter_force,
            time_step=min(COBRA_TIMESTEP, BIOSCRAPE_TIMESTEP), #/2
        )
        lattice_config['multibody']['_parallel'] = parallel
        lattice_composer = Lattice(lattice_config)

        # merge bioscrapeCOBRA composite with lattice
        biocobra_composite.merge(composite=lattice_composer.generate())
        initial_composite.merge(composite=lattice_composer.generate())

        # get initial state from the composite and merge declared initial
        state = biocobra_composite.initial_state()
        initial_state_full = deep_merge(state, agents_initial)

    elif division:
        # make n_agents dividing agents, without an explicit environment

        # division requires the agent to be embedded in a hierarchy
        # make the bioscrapeCOBRA composite under the path ('agents', agent_id)
        agents_initial = {'agents': {}}
        for n in range(n_agents):
            agent_id = str(n)
            config = {'agent_id': agent_id}
            agent = biocobra_composer.generate(config=config)
            biocobra_composite.merge(
                composite=agent,
                path=('agents', agent_id))

            # initial state for agent
            agents_initial['agents'][agent_id] = get_initial_state(initial_agent_states, n=n)
            agents_initial['agents'][agent_id] = deep_merge(
                agents_initial['agents'][agent_id], agent_state)

        # create a second initial composite for plotting
        initial_composite = biocobra_composer.generate(
            path=('agents', '0'),
            config={'agent_id': '0'})

        # get initial state from the composite and merge declared initial
        state = biocobra_composite.initial_state()
        initial_state_full = deep_merge(state, agents_initial)

        # simple, 1D field
        initial_state_full['fields'] = {
                GLUCOSE_EXTERNAL: initial_glucose,
                LACTOSE_EXTERNAL: initial_lactose}

    else:
        # single agent without division

        # make the composite
        agent = biocobra_composer.generate()
        biocobra_composite.merge(
            composite=agent,
            path=())

        # create a second initial composite for plotting
        initial_composite = biocobra_composer.generate()

        # get initial state from the composite and merge declared initial
        agents_initial = biocobra_composite.initial_state()
        initial_declared = get_initial_state(initial_agent_states)
        agents_initial = deep_merge(agents_initial, initial_declared)
        initial_state_full = deep_merge(agents_initial, agent_state)

        # simple, 1D field
        initial_state_full['fields'] = {
                GLUCOSE_EXTERNAL: initial_glucose,
                LACTOSE_EXTERNAL: initial_lactose}

    # make the experiment
    experiment_id = (f"{'stochastic' if stochastic else 'deterministic'}"
                     f"{'_division' if division else ''}"
                     f"{'_spatial' if spatial else ''}"
                     f"_{timestamp()}")
    experiment_config = {
        'processes': biocobra_composite.processes,
        'topology': biocobra_composite.topology,
        'initial_state': initial_state_full,
        'display_info': False,
        'experiment_id': experiment_id,
        'emit_step': max(BIOSCRAPE_TIMESTEP, COBRA_TIMESTEP),
        'emitter': {'type': emitter}}
    print(f'Initializing experiment {experiment_id}')
    biocobra_experiment = Engine(**experiment_config)

    # run the experiment
    clock_start = clock.time()
    if division:  # terminate upon reaching total_time or halt_threshold
        sim_step = max(BIOSCRAPE_TIMESTEP, COBRA_TIMESTEP) * 10
        for _ in tqdm(range(0, total_time, sim_step)):
            n_agents = len(biocobra_experiment.state.get_value()['agents'])
            if n_agents < halt_threshold:
                biocobra_experiment.update(sim_step)
    else:
        biocobra_experiment.update(total_time)

    # print runtime and finalize
    clock_finish = clock.time() - clock_start
    print(f'Completed in {clock_finish:.2f} seconds')
    biocobra_experiment.end()

    # retrieve the data
    if output_type == 'timeseries':
        return biocobra_experiment.emitter.get_timeseries(), initial_composite
    if output_type == 'unitless':
        return biocobra_experiment.emitter.get_data_unitless(), initial_composite
    return biocobra_experiment, initial_composite


# plot topology
def plot_full_topology(out_dir='out'):

    dirname = os.path.dirname(__file__)
    sbml_file = os.path.join(dirname, SBML_FILE_STOCHASTIC)

    # get the stochastic bioscrape_cobra composite
    output, composite = simulate_bioscrape_cobra(
        stochastic=True,
        division=True,
        spatial=True,
        total_time=0,
        sbml_file=sbml_file,
        output_type='unitless')

    plot_topology(
        composite,
        settings=config_embedded_bioscrape_cobra_topology,
        out_dir=out_dir,
        filename='bioscrape_cobra_stochastic_lattice_topology.pdf')


# tests
def test_deterministic(
        emitter=None,
        sbml_deterministic=None,
        total_time=200,
        out_dir=None,
):
    if not sbml_deterministic:
        dirname = os.path.dirname(__file__)
        sbml_deterministic = os.path.join(dirname, SBML_FILE_DETERMINISTIC)

    output, comp0 = simulate_bioscrape_cobra(
        initial_glucose=1e1,
        initial_lactose=2e1,
        external_volume=1e-12,
        total_time=total_time,
        emitter=emitter,
        sbml_file=sbml_deterministic,
        output_type='timeseries')

    if out_dir:
        plot_single(
            output,
            variables=plot_variables_list_deterministic,
            out_dir=os.path.join(out_dir, 'deterministic'),
            filename='variables')


@pytest.mark.slow
def test_stochastic(
        emitter=None,
        sbml_stochastic=None,
        total_time=500,
        out_dir=None,
):
    if not sbml_stochastic:
        dirname = os.path.dirname(__file__)
        sbml_stochastic = os.path.join(dirname, SBML_FILE_STOCHASTIC)

    initial_agent_state = {
        'rates': {
            'k_leak': 0.6,
        },
        'species': {
            'monomer_betaGal': 100,
            'protein_betaGal': 100,
            'protein_Lactose_Permease': 100}}

    output, comp0 = simulate_bioscrape_cobra(
        stochastic=True,
        initial_glucose=1e1,
        initial_lactose=2e1,
        initial_agent_states=initial_agent_state,
        external_volume=1e-14,
        total_time=total_time,
        emitter=emitter,
        sbml_file=sbml_stochastic,
        output_type='timeseries')

    if out_dir:
        plot_single(
            output,
            variables=plot_variables_list_stochastic,
            out_dir=os.path.join(out_dir, 'stochastic'),
            filename='variables')


def test_deterministic_divide(
        emitter=None,
        sbml_deterministic=None,
        total_time=200,
        out_dir=None,
):
    if not sbml_deterministic:
        dirname = os.path.dirname(__file__)
        sbml_deterministic = os.path.join(dirname, SBML_FILE_DETERMINISTIC)

    initial_agent_states = [
        {'species': {
            'monomer_betaGal': 0.0,
            'protein_betaGal': 0.0,
            'protein_Lactose_Permease': 0.0}},
        {'species': {
            'monomer_betaGal': 0.0,
            'protein_betaGal': 0.1,
            'protein_Lactose_Permease': 0.1}}]

    output, comp0 = simulate_bioscrape_cobra(
        n_agents=2,
        initial_agent_states=initial_agent_states,
        division=True,
        initial_glucose=1e1,
        initial_lactose=2e1,
        external_volume=1e-12,
        total_time=total_time,
        emitter=emitter,
        sbml_file=sbml_deterministic,
        output_type='unitless')

    if out_dir:
        var_list = copy.deepcopy(plot_variables_list_deterministic)
        plot_multigen(
            output,
            variables=var_list,
            out_dir=os.path.join(out_dir, 'deterministic_divide'),
            filename='division_multigen')


@pytest.mark.slow
def test_stochastic_divide(
        emitter=None,
        sbml_stochastic=None,
        total_time=200,
        out_dir=None,
):
    if not sbml_stochastic:
        dirname = os.path.dirname(__file__)
        sbml_stochastic = os.path.join(dirname, SBML_FILE_STOCHASTIC)

    initial_agent_states = [{
            'rates': {
                'LacPermease_vmax': 3580.0,  # 35.8
                'k_leak': 0.5,
            },
        }]

    output, comp0 = simulate_bioscrape_cobra(
        n_agents=1,
        stochastic=True,
        division=True,
        initial_glucose=1e1,  # mM
        initial_lactose=2e1,  # mM
        initial_agent_states=initial_agent_states,
        total_time=total_time,
        external_volume=5e-14,
        emitter=emitter,
        sbml_file=sbml_stochastic,
        output_type='unitless')

    if out_dir:
        # plot
        var_list = copy.deepcopy(plot_variables_list_stochastic)
        var_list.extend([
            ('boundary', 'mass'),
            ('boundary', 'volume')])
        plot_multigen(
            output,
            variables=var_list,
            out_dir=os.path.join(out_dir, 'stochastic_divide'),
            filename='division_multigen')


def test_deterministic_spatial(
        emitter=None,
        sbml_deterministic=None,
        total_time=200,
        out_dir=None,
):
    if not sbml_deterministic:
        dirname = os.path.dirname(__file__)
        sbml_deterministic = os.path.join(dirname, SBML_FILE_DETERMINISTIC)

    bounds = [30, 30]
    n_bins = [30, 30]
    depth = 1

    initial_agent_states = [
        {'rates': {
            'k_leak': 0.5,
        }}
    ]

    output, comp0 = simulate_bioscrape_cobra(
        initial_agent_states=initial_agent_states,
        division=True,
        spatial=True,
        initial_glucose=1e1,
        initial_lactose=2e1,
        bounds=bounds,
        n_bins=n_bins,
        depth=depth,
        total_time=total_time,
        emitter=emitter,
        sbml_file=sbml_deterministic,
        output_type='unitless')

    if out_dir:
        # plot
        deterministic_spatial_out_dir = os.path.join(out_dir, 'deterministic_spatial')
        plot_multigen(
            output,
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_multigen')

        plot_fields_snapshots(
            output,
            bounds=bounds,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_snapshots')

        plot_fields_tags(
            output,
            bounds=bounds,
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_tags')


@pytest.mark.slow
def test_stochastic_spatial(
    emitter=None,
    parallel=False,
    sbml_stochastic=None,
    total_time=200,
    out_dir=None,
):
    if not sbml_stochastic:
        dirname = os.path.dirname(__file__)
        sbml_stochastic = os.path.join(dirname, SBML_FILE_STOCHASTIC)

    bounds = [30, 30]
    n_bins = [30, 30]
    depth = 0.5

    initial_agent_states = [
        {'rates': {
            'k_leak': 0.005  # less leak -> less spontanteous expression
        }}
    ]

    output, comp0 = simulate_bioscrape_cobra(
        initial_agent_states=initial_agent_states,
        stochastic=True,
        division=True,
        spatial=True,
        initial_glucose=1e1,
        initial_lactose=5e1,
        depth=depth,
        diffusion_rate=2e-2,
        jitter_force=1e-5,
        bounds=bounds,
        n_bins=n_bins,
        halt_threshold=200,
        total_time=total_time,
        emitter=emitter,
        sbml_file=sbml_stochastic,
        parallel=parallel,
        output_type='unitless')

    if out_dir:
        # plot
        stochastic_spatial_out_dir = os.path.join(out_dir, 'stochastic_spatial')
        plot_multigen(
            output,
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_multigen')

        plot_fields_snapshots(
            output,
            bounds=bounds,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_snapshots')

        plot_fields_tags(
            output,
            bounds=bounds,
            convert_to_concs=False,
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_tags')


def test_bioscrape_alone(
        sbml_deterministic=None,
        sbml_stochastic=None,
        out_dir=None,
):
    if not sbml_deterministic:
        dirname = os.path.dirname(__file__)
        sbml_deterministic = os.path.join(dirname, SBML_FILE_DETERMINISTIC)
    if not sbml_stochastic:
        dirname = os.path.dirname(__file__)
        sbml_stochastic = os.path.join(dirname, SBML_FILE_STOCHASTIC)


    # Simulate the Lac Operon CRN Deterministically
    bioscrape_timeseries_det, bioscrape_composite_det = simulate_bioscrape(
        total_time=20000,
        initial_glucose=10,
        initial_lactose=20,
        sbml_file=sbml_deterministic)
    # Simulate the Lac Operon CRN Stochastically
    bioscrape_timeseries_stoch, bioscrape_composite_stoch = simulate_bioscrape(
        total_time=20000,
        initial_glucose=10 * 6.22 * 10 ** 6,
        initial_lactose=20 * 6.22 * 10 ** 6,
        stochastic=True,
        sbml_file=sbml_stochastic)

    if out_dir:
        # Plot the CRN Trajectories
        species_to_plot_det = [
            {'variable': ('species', 'Glucose_external'), 'display': 'external glucose (mM)'},
            {'variable': ('species', 'Lactose_external'), 'display': 'external lactose (mM)'},
            {'variable': ('species', 'rna_M'), 'display': 'lac operon RNA (mM)'},
            {'variable': ('species', 'protein_betaGal'), 'display': r'$\beta$-galactosidase (mM)'},
            {'variable': ('species', 'protein_Lactose_Permease'), 'display': 'lactose permease (mM)'}
        ]
        species_to_plot_stoch = [
            {'variable': ('species', 'Glucose_external'), 'display': 'external glucose  (counts)'},
            {'variable': ('species', 'Lactose_external'), 'display': 'external lactose (counts)'},
            {'variable': ('species', 'rna_M'), 'display': 'lac operon RNA (counts)'},
            {'variable': ('species', 'protein_betaGal'), 'display': r'$\beta$-galactosidase (counts)'},
            {'variable': ('species', 'protein_Lactose_Permease'), 'display': 'lactose permease (counts)'}
        ]

        fig_timeseries = plot_single(
            bioscrape_timeseries_det,
            variables=species_to_plot_det)
        save_fig_to_dir(
            fig_timeseries,
            filename='bioscrape_deterministic_output.pdf',
            out_dir=out_dir)

        fig_timeseries = plot_single(
            bioscrape_timeseries_stoch,
            variables=species_to_plot_stoch)
        save_fig_to_dir(
            fig_timeseries,
            filename='bioscrape_stochastic_output.pdf',
            out_dir=out_dir)


def main():
    out_dir = os.path.join('out', 'bioscrape_cobra')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='bioscrape_cobra')
    parser.add_argument('-database', '-d', action='store_true', default=False, help='emit to database')
    parser.add_argument('-parallel', '-p', action='store_true', default=False, help='run parallel processes')
    parser.add_argument('--deterministic', '-1', action='store_true', default=False)
    parser.add_argument('--stochastic', '-2', action='store_true', default=False)
    parser.add_argument('--deterministic_divide', '-3', action='store_true', default=False)
    parser.add_argument('--stochastic_divide', '-4', action='store_true', default=False)
    parser.add_argument('--deterministic_spatial', '-5', action='store_true', default=False)
    parser.add_argument('--stochastic_spatial', '-6', action='store_true', default=False)
    parser.add_argument('--topology', '-t', action='store_true', default=False)
    parser.add_argument('--bioscrape-alone', '-b', action='store_true', default=False)

    args = parser.parse_args()

    # emitter type
    emitter = 'database' if args.database else 'timeseries'
    parallel = True if args.parallel else False

    if args.deterministic:
        test_deterministic(
            emitter=emitter,
            sbml_deterministic=sbml_deterministic_file,
            total_time=12000,
            out_dir=out_dir
        )

    if args.stochastic:
        test_stochastic(
            emitter=emitter,
            sbml_stochastic=sbml_stochastic_file,
            total_time=4000,
            out_dir=out_dir,
        )

    if args.deterministic_divide:
        test_deterministic_divide(
            emitter=emitter,
            sbml_deterministic=sbml_deterministic_file,
            total_time=3000,
            out_dir=out_dir
        )

    if args.stochastic_divide:
        test_stochastic_divide(
            emitter=emitter,
            sbml_stochastic=sbml_stochastic_file,
            total_time=4000,
            out_dir=out_dir,
        )

    if args.deterministic_spatial:
        test_deterministic_spatial(
            emitter=emitter,
            sbml_deterministic=sbml_deterministic_file,
            total_time=7200,
            out_dir=out_dir,
        )

    if args.stochastic_spatial:
        test_stochastic_spatial(
            emitter=emitter,
            parallel=parallel,
            sbml_stochastic=sbml_stochastic_file,
            total_time=50400,
            out_dir=out_dir,
        )

    if args.topology:
        plot_full_topology(out_dir=out_dir)

    if args.bioscrape_alone:
        test_bioscrape_alone()


if __name__ == '__main__':
    main()

