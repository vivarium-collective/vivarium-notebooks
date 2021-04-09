"""
==============================================
Simulation helper functions for BioscrapeCOBRA
==============================================
"""
import os
import argparse
import copy
import random
import time as clock
from tqdm import tqdm

# vivarium imports
from vivarium.core.experiment import Experiment, timestamp
from vivarium.library.units import units
from vivarium.library.dict_utils import deep_merge
from vivarium.core.composition import simulate_composer, simulate_composite
from vivarium.core.process import Composite
from vivarium.library.units import remove_units
from vivarium.core.process import deserialize_value

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import (
    Lattice, make_lattice_config)
from vivarium_multibody.composites.grow_divide import GrowDivide
from vivarium_multibody.processes.diffusion_field import DiffusionField

# cobra imports
from vivarium_cobra.processes.dynamic_fba import (
    DynamicFBA, 
    get_iAF1260b_config, 
)
from vivarium_cobra.composites.cobra_composite import CobraComposite
from vivarium_cobra.library.lattice_utils import get_bin_volume

# Bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# local import
from bioscrape_cobra.bioscrape_cobra_stochastic import BioscrapeCOBRAstochastic, SBML_FILE_STOCHASTIC
from bioscrape_cobra.bioscrape_cobra_stochastic import GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL
from bioscrape_cobra.bioscrape_cobra_deterministic import BioscrapeCOBRAdeterministic, SBML_FILE_DETERMINISTIC

# plotting
from bioscrape_cobra.plot import (
    plot_multigen, plot_single, plot_fields_tags, plot_fields_snapshots)

# default variables, which can be varied by simulate_bioscrape_cobra
DEFAULT_DIVIDE_THRESHOLD = 2000 * units.fg
INITIAL_GLC = 1e2  # mmolar
INITIAL_LAC = 1e2  # mmolar
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 2

# fixed global variables
COBRA_TIMESTEP = 50
BIOSCRAPE_TIMESTEP = 10


# Simulates Cobra Model on its Own
def simulate_cobra(
    total_time=100,
    initial_state=None,
):
    # get the configuration for the iAF1260b BiGG model
    iAF1260b_config = get_iAF1260b_config()
    iAF1260b_config.update({'time_step': COBRA_TIMESTEP})

    iAF1260b_config['name'] = 'cobra' #Rename the process
    dynamic_fba = DynamicFBA(iAF1260b_config)

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
        initial_volume=1.0
    ):
    
    #create configs
    if not stochastic:
        bioscrape_config = {
            'sbml_file': 'LacOperon_deterministic.xml',
            'stochastic': False,
            'initial_volume': initial_volume,
            'internal_dt': BIOSCRAPE_TIMESTEP/100,
            'time_step':BIOSCRAPE_TIMESTEP}
    else:
        bioscrape_config = {
            'sbml_file': 'LacOperon_stochastic.xml',
            'stochastic': True,
            'safe_mode': False,
            'initial_volume': initial_volume,
            'internal_dt': BIOSCRAPE_TIMESTEP/100,
            'time_step':BIOSCRAPE_TIMESTEP
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


# Simulate a System of Cells that grow and divde in a well mixed spatial environment
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

# Simulate a System of Cells that grow and divde in a well mixed spatial environment
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

    if initial_state is None:
        initial_state = diffusion_process.initial_state

    diffusion_sim_settings = {
        'total_time': total_time,
        'return_raw_data': True,
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
        bounds=[20, 20],
        n_bins=[10, 10],
        depth=10,
        diffusion_rate=1e-1,
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
    """ Simulation function for BioscrapeCOBRA

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

    # get the bin volume based upon the lattice
    bin_volume = (external_volume or get_bin_volume(n_bins, bounds, depth)) * units.L
    agent_state = {
        # set initial flux values based on COBRA defaults.
        'flux_bounds': {
            'EX_lac__D_e': 0.0,
            'EX_glc__D_e': 0.099195},
        # field_counts_deriver needs the bin volume
        'boundary': {
            'bin_volume': bin_volume,
            **({'external': {
                GLUCOSE_EXTERNAL: initial_glucose,
                LACTOSE_EXTERNAL: initial_lactose}} if spatial else {}),
            # 'external': {
            #     GLUCOSE_EXTERNAL: initial_glucose,
            #     LACTOSE_EXTERNAL: initial_lactose}
        }}

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

    # make the composite
    if spatial:
        # make n_agents dividing agents and combine them with a Lattice environment

        # make a bioscrapeCOBRA composite
        biocobra_composite = Composite({})
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

        # get initial state from the composite and merge declared initial
        state = biocobra_composite.initial_state()
        initial_state_full = deep_merge(state, agents_initial)

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
            time_step=min(COBRA_TIMESTEP, BIOSCRAPE_TIMESTEP))
        lattice_config['multibody']['_parallel'] = parallel

        lattice_composer = Lattice(lattice_config)
        lattice_composite = lattice_composer.generate()

        # merge bioscrapeCOBRA composite with lattice
        biocobra_composite.merge(composite=lattice_composite)
        initial_composite.merge(composite=lattice_composer.generate())

    elif division:
        # make n_agents dividing agents, without an explicit environment

        # division requires the agent to be embedded in a hierarchy
        # make the bioscrapeCOBRA composite under the path ('agents', agent_id)
        biocobra_composite = Composite({})
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
        biocobra_composite = biocobra_composer.generate()

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
    biocobra_experiment = Experiment(experiment_config)

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


# plotting
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

def main():
    out_dir = os.path.join(
        'out', 'bioscrape_cobra')
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
    args = parser.parse_args()

    # emitter type
    emitter = 'database' if args.database else 'timeseries'
    parallel = True if args.parallel else False

    # set sbml file path
    dirname = os.path.dirname(__file__)
    sbml_deterministic = os.path.join(dirname, SBML_FILE_DETERMINISTIC)
    sbml_stochastic = os.path.join(dirname, SBML_FILE_STOCHASTIC)

    if args.deterministic:
        output, comp0 = simulate_bioscrape_cobra(
            initial_glucose=1e0,
            initial_lactose=1e0,
            external_volume=1e-13,
            total_time=3000,
            emitter=emitter,
            sbml_file=sbml_deterministic,
            output_type='timeseries')

        plot_single(
            output,
            variables=plot_variables_list_deterministic,
            out_dir=os.path.join(out_dir, 'deterministic'),
            filename='variables')

    if args.stochastic:
        output, comp0 = simulate_bioscrape_cobra(
            stochastic=True,
            initial_glucose=1e0,
            initial_lactose=1e0,
            external_volume=1e-14,
            total_time=3000,
            emitter=emitter,
            sbml_file=sbml_stochastic,
            output_type='timeseries')

        plot_single(
            output,
            variables=plot_variables_list_stochastic,
            out_dir=os.path.join(out_dir, 'stochastic'),
            filename='variables')

    if args.deterministic_divide:
        initial_agent_states = [{
            'species': {
                'monomer_betaGal': 0.0,
                'protein_betaGal': 0.0,
                'protein_Lactose_Permease': 0.0},
            },
            {
                'species': {
                    'monomer_betaGal': 0.0,
                    'protein_betaGal': 0.1,
                    'protein_Lactose_Permease': 0.1},
            },
        ]

        output, comp0 = simulate_bioscrape_cobra(
            n_agents=2,
            initial_agent_states=initial_agent_states,
            division=True,
            initial_glucose=1e1,  # mM
            initial_lactose=1e1,  # mM
            external_volume=1e-12,
            total_time=4000,
            emitter=emitter,
            sbml_file=sbml_deterministic,
            output_type='unitless')

        var_list = copy.deepcopy(plot_variables_list_deterministic)
        plot_multigen(
            output,
            variables=var_list,
            out_dir=os.path.join(out_dir, 'deterministic_divide'),
            filename='division_multigen')

    if args.stochastic_divide:
        initial_agent_states = [{
            'species': {
                'monomer_betaGal': 0,
                'protein_betaGal': 0,
                'protein_Lactose_Permease': 0}
        },
        {
            'species': {
                'monomer_betaGal': 100,
                'protein_betaGal': 100,
                'protein_Lactose_Permease': 100}
        }]

        output, comp0 = simulate_bioscrape_cobra(
            n_agents=2,
            stochastic=True,
            division=True,
            initial_glucose=1e0,  # mM
            initial_lactose=1e1,  # mM
            initial_agent_states=initial_agent_states,
            total_time=4000,
            external_volume=1e-12,
            emitter=emitter,
            sbml_file=sbml_stochastic,
            output_type='unitless')

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

    if args.deterministic_spatial:

        output, comp0 = simulate_bioscrape_cobra(
            division=True,
            spatial=True,
            initial_glucose=1e1,
            initial_lactose=1e1,
            total_time=12000,
            emitter=emitter,
            sbml_file=sbml_deterministic,
            output_type='unitless')

        deterministic_spatial_out_dir = os.path.join(out_dir, 'deterministic_spatial')
        plot_multigen(
            output,
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_multigen')

        plot_fields_snapshots(
            output,
            bounds=BOUNDS,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_snapshots')

        plot_fields_tags(
                        output,
            bounds=BOUNDS,
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_tags')

    if args.stochastic_spatial:
        bounds = [25, 25]
        n_bins = [25, 25]

        output, comp0 = simulate_bioscrape_cobra(
            stochastic=True,
            division=True,
            spatial=True,
            initial_glucose=1e1,
            initial_lactose=1e5,
            depth=1,
            diffusion_rate=1e-1,
            initial_state=None,
            bounds=bounds,
            n_bins=n_bins,
            halt_threshold=200,
            total_time=60000,
            emitter=emitter,
            sbml_file=sbml_stochastic,
            parallel=parallel,
            output_type='unitless')

        stochastic_spatial_out_dir = os.path.join(out_dir, 'stochastic_spatial')
        plot_multigen(
            output,
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_multigen')

        plot_fields_snapshots(
            output,
            bounds=bounds,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            n_snapshots=5,
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_snapshots')

        plot_fields_tags(
            output,
            bounds=bounds,
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            n_snapshots=5,
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_tags')


if __name__ == '__main__':
    main()
