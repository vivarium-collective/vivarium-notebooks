"""
==============================================
Simulation helper functions for BioscrapeCOBRA
==============================================
"""
import os
import argparse

# vivarium imports
from vivarium.core.experiment import Experiment
from vivarium.library.units import units
from vivarium.library.dict_utils import deep_merge
from vivarium.core.composition import simulate_composer, composer_in_experiment, simulate_composite
from vivarium.core.process import Composite
from vivarium.library.units import remove_units
from vivarium.core.process import deserialize_value

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import (
    Lattice, make_lattice_config)
from vivarium_multibody.processes.multibody_physics import test_growth_division
from vivarium_multibody.processes.multibody_physics import agent_body_config, volume_from_length
from vivarium_multibody.composites.grow_divide import GrowDivide
from vivarium_multibody.plots.snapshots import (
    plot_snapshots,
    format_snapshot_data,
    make_snapshots_figure,
    get_field_range,
    get_agent_colors,
)
from vivarium_multibody.processes.diffusion_field import DiffusionField, get_gaussian_config

#cobra imports
from vivarium_cobra.processes.dynamic_fba import (
    DynamicFBA, 
    get_iAF1260b_config, 
    print_growth
)
from vivarium_cobra.composites.cobra_composite import CobraComposite
from vivarium_cobra.library.lattice_utils import get_bin_volume

#Bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# local import
from bioscrape_cobra.bioscrape_cobra_stochastic import (
    BioscrapeCOBRAstochastic, GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL)
from bioscrape_cobra.bioscrape_cobra_deterministic import BioscrapeCOBRAdeterministic

# plotting
from bioscrape_cobra.plot import (plot_multigen, plot_single, plot_fields_tags, plot_fields_snapshots)

# default variables, which can be varied by simulate_bioscrape_cobra
DEFAULT_EXTERNAL_VOLUME = 1e-13 * units.L
DEFAULT_DIVIDE_THRESHOLD = 2000 * units.fg
INITIAL_GLC = 10e1  # mmolar
INITIAL_LAC = 10e1  # mmolar
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 20

# fixed global variables
COBRA_TIMESTEP = 50
BIOSCRAPE_TIMESTEP = 10

# divide config
INITIAL_AGENT_ID = '1'
divide_config = {
    'divide_on': True,
    'agent_id': INITIAL_AGENT_ID,  # TODO -- this should be configured below in case it is overwritten
    'agents_path': ('..', '..', 'agents',),
    'fields_path': ('..', '..', 'fields',),
    'dimensions_path': ('..', '..', 'dimensions',),
    'local_fields': {}}

# spatial config
spatial_config = dict(divide_config)
spatial_config['fields_on'] = True

#Simulates Cobra Model on its Own
def simulate_cobra(
    total_time=100,
    initial_state = None):
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
        'total_time': total_time
        }

    cobra_timeseries = simulate_composer(dynamic_fba, cobra_sim_settings)

    return cobra_timeseries, dynamic_fba

#simulates Cobra Model in a Composite with some Derivers
def simulate_cobra_composite(
    total_time=100,
    initial_state = None):

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
        'total_time': total_time
        }

    cobra_timeseries = simulate_composer(cobra_composite, cobra_sim_settings)

    return cobra_timeseries, cobra_composite

def simulate_bioscrape(
        stochastic=False,
        initial_glucose=None,
        initial_lactose=None,
        initial_state=None,
        total_time=100,
        initial_volume = 1.0
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
        initial_state["species"]["Glucose_external"] = initial_glucose

    if initial_lactose is not None:
        initial_state["species"]["Lactose_external"] = initial_lactose

    # run simulation
    bioscrape_sim_settings = {
        'initial_state': initial_state,
        'total_time': total_time
        }

    bioscrape_timeseries = simulate_composer(bioscrape_process, bioscrape_sim_settings)

    return bioscrape_timeseries, bioscrape_composite


#Simulate a System of Cells that grow and divde in a well mixed spatial environment
def simulate_grow_divide(total_time = 100, growth_rate = .03, initial_state = None, growth_noise = 10**-6):
    # config
    growth_config = {'default_growth_rate': growth_rate, "default_growth_noise": growth_noise, 'agents_path': ('agents',)}
    grow_divide_composer = GrowDivide({'agent_id': "0", 'growth' : growth_config})

    if initial_state is None:
        initial_state = {
        'agents': {
            '0': {
                'global': {
                    'mass': 1000 * units.femtogram}
            }}}

    grow_divide_sim_settings = {
        'total_time': total_time,
        "initial_state": initial_state,
        "outer_path":('agents', '0'),
        'return_raw_data': True,
    }

    grow_divide_composite = grow_divide_composer.generate(path=('agents', "0"))

    grow_divide_data = simulate_composite(grow_divide_composite, grow_divide_sim_settings)
    grow_divide_data = deserialize_value(grow_divide_data)
    grow_divide_data = remove_units(grow_divide_data)

    #returns the data, the initial composite, and the final composite
    return grow_divide_data, grow_divide_composer.generate(path=('agents', "0")), grow_divide_composite

#Simulate a System of Cells that grow and divde in a well mixed spatial environment
def simulate_diffusion(total_time = 100, diffusion_rate = .001, initial_state = {}, bins = [10, 10], bounds = [10, 10]):
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


def get_lattice_grow_divide_composite(diffusion_rate = .001, initial_concentration = {}, bins = [10, 10], bounds = [10, 10], growth_rate = .03, growth_noise = 10**-6, depth = 10):
    lattice_config = make_lattice_config(
            bounds=bounds,
            n_bins=bins,
            concentrations={'glc':initial_concentration},
            diffusion=diffusion_rate,
            depth = 10)

    lattice_composer = Lattice(lattice_config)
    lattice_composite = lattice_composer.generate()


    growth_config = {'default_growth_rate': growth_rate, "default_growth_noise": growth_noise}
    grow_divide_composer = GrowDivide({'agent_id': "0", 'growth' : growth_config})

    agent_id = "0"
    lattice_grow_divide_composite = grow_divide_composer.generate(path=('agents', agent_id))
    lattice_grow_divide_composite.merge(composite=lattice_composite)

    return lattice_grow_divide_composite

def simulate_grow_divide_lattice(lattice_grow_divide_composite, total_time = 100,  initial_state = None):

    agent_id = "0"
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
        "initial_state" : initial_state,
        'return_raw_data': True,
    }

    lattice_grow_divide_data = simulate_composite(lattice_grow_divide_composite, sim_settings)

    return lattice_grow_divide_data


# helper functions
def get_bioscrape_cobra_config(
        spatial=False,
        division=False,
        divide_threshold=DEFAULT_DIVIDE_THRESHOLD,
        external_volume=DEFAULT_EXTERNAL_VOLUME,
        agent_id=INITIAL_AGENT_ID
):
    """ create a generic config dict for bioscrape_cobra composers """
    agent_id = agent_id

    if spatial:
        config = {
            'divide_on': True,
            'fields_on': True,
            '_parallel': True,
            'agent_id': agent_id,
            'agents_path': ('..', '..', 'agents',),
            'fields_path': ('..', '..', 'fields',),
            'dimensions_path': ('..', '..', 'dimensions',),
            'local_fields': {},
            'divide_condition': {
                'threshold': divide_threshold}}
    elif division:
        config = {
            'divide_on': True,
            '_parallel': True,
            'agent_id': agent_id,
            'agents_path': ('..', '..', 'agents',),
            'fields_path': ('..', '..', 'fields',),
            'dimensions_path': ('..', '..', 'dimensions',),
            'local_fields': {},
            'divide_condition': {
                'threshold': divide_threshold}}
    else:
        config = {
            'local_fields': {
                'bin_volume': external_volume},
            '_parallel': False}

    return config



def simulate_bioscrape_cobra(
        division=False,
        stochastic=False,
        spatial=False,
        initial_glucose=1e1,
        initial_lactose=1e1,
        initial_state=None,
        bounds=[20, 20],
        n_bins=[10, 10],
        depth=10,
        diffusion_rate=1e-1,
        divide_threshold=2000 * units.fg,
        external_volume=1e-13 * units.L,
        agent_id='1',
        halt_threshold=32,
        total_time=100,
        output_type=None,
):
    """ Simulation function for BioscrapeCOBRA """
    initial_state = initial_state or {}

    # make the BioscrapeCOBRA config
    biocobra_config = get_bioscrape_cobra_config(
        spatial=spatial,
        division=division,
        divide_threshold=divide_threshold,
        external_volume=external_volume,
        agent_id=agent_id)

    # get the BioscrapeCOBRA composer -- either stochastic or deterministic
    if stochastic:
        biocobra_composer = BioscrapeCOBRAstochastic(biocobra_config)
    else:
        biocobra_composer = BioscrapeCOBRAdeterministic(biocobra_config)

    # make the composite
    if spatial:
        # make a bioscrapeCOBRA composite
        biocobra_composite = biocobra_composer.generate(
            path=('agents', agent_id))

        # get initial state from the composite

        #set the bin volume based upon the lattice
        bin_volume = get_bin_volume(n_bins, bounds, depth)
        bin_volume_config = config = {'local_fields': {'bin_volume': bin_volume}}
        state = biocobra_composite.initial_state(bin_volume_config)
        state['agents'][agent_id]['boundary']['external'] = {
            GLUCOSE_EXTERNAL: initial_glucose,
            LACTOSE_EXTERNAL: initial_lactose}

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
            time_step=COBRA_TIMESTEP)
        lattice_config['multibody']['_parallel'] = True
        lattice_config['multibody']['timestep'] = BIOSCRAPE_TIMESTEP

        lattice_composer = Lattice(lattice_config)
        lattice_composite = lattice_composer.generate()

        # merge bioscrapeCOBRA composite with lattice
        biocobra_composite.merge(composite=lattice_composite)

    elif division:
        # division requires the agent to be embedded in a hierarchy
        # make the bioscrapeCOBRA composite under the path ('agents', agent_id)
        biocobra_composite = biocobra_composer.generate(
            path=('agents', agent_id))

        # get initial state from the composite
        state = biocobra_composite.initial_state()
        state['agents'][agent_id]['boundary']['external'] = {
            GLUCOSE_EXTERNAL: initial_glucose,
            LACTOSE_EXTERNAL: initial_lactose}

    else:
        # make the composite
        biocobra_composite = biocobra_composer.generate()

        # get initial state from the composite
        state = biocobra_composite.initial_state()
        state['boundary']['external'] = {
            GLUCOSE_EXTERNAL: initial_glucose,
            LACTOSE_EXTERNAL: initial_lactose}

    # update initial state with any value from function assignment
    initial_state = deep_merge(state, initial_state)

    # make the experiment
    experiment_config = {
        'processes': biocobra_composite.processes,
        'topology': biocobra_composite.topology,
        'initial_state': initial_state,
        # 'display_info': False,
        'experiment_id': f"{'stochastic' if stochastic else 'deterministic'}"
                         f"{'_division' if division else ''}"
                         f"{'_spatial' if spatial else ''}"}
    biocobra_experiment = Experiment(experiment_config)

    # run the experiment
    if division:
        # terminate upon reaching total_time or halt_threshold
        time = 0
        sim_step = 1000
        n_agents = len(biocobra_experiment.state.get_value()['agents'])
        while n_agents < halt_threshold and time <= total_time:
            biocobra_experiment.update(sim_step)
            time += sim_step
            n_agents = len(biocobra_experiment.state.get_value()['agents'])
    else:
        biocobra_experiment.update(total_time)

    # retrieve the data
    if output_type == 'timeseries':
        return biocobra_experiment.emitter.get_timeseries()
    if output_type == 'unitless':
        return biocobra_experiment.emitter.get_data_unitless()
    return biocobra_experiment


# plotting
plot_variables_list = [
    ('species', 'rna_M'),
    ('species', 'protein_betaGal'),
    ('species', 'protein_Lactose_Permease'),
    ('flux_bounds', 'EX_glc__D_e'),
    ('flux_bounds', 'EX_lac__D_e'),
    ('boundary', ('mass', 'femtogram')),
    ('boundary', ('volume', 'femtoliter'))]

plot_variables_list_deterministic = [
    ('boundary', 'external', GLUCOSE_EXTERNAL),
    ('boundary', 'external', LACTOSE_EXTERNAL),
    ('rates', 'k_dilution__',)]
plot_variables_list_deterministic.extend(plot_variables_list)

plot_variables_list_stochastic = [
    ('species', GLUCOSE_EXTERNAL),
    ('species', LACTOSE_EXTERNAL)]
plot_variables_list_stochastic.extend(plot_variables_list)

def main():
    out_dir = os.path.join(
        'out', 'bioscrape_cobra')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='bioscrape_cobra')
    parser.add_argument('--deterministic', '-1', action='store_true', default=False)
    parser.add_argument('--stochastic', '-2', action='store_true', default=False)
    parser.add_argument('--deterministic_divide', '-3', action='store_true', default=False)
    parser.add_argument('--stochastic_divide', '-4', action='store_true', default=False)
    parser.add_argument('--deterministic_spatial', '-5', action='store_true', default=False)
    parser.add_argument('--stochastic_spatial', '-6', action='store_true', default=False)
    args = parser.parse_args()

    if args.deterministic:
        output = simulate_bioscrape_cobra(
            total_time=2000,
            output_type='timeseries')

        plot_single(
            output,
            variables=plot_variables_list_deterministic,
            out_dir=os.path.join(out_dir, 'deterministic'),
            filename='variables')

    if args.stochastic:
        initial_state = {
            'species': {
                'protein_Lactose_Permease': 10}}

        output = simulate_bioscrape_cobra(
            stochastic=True,
            initial_state=initial_state,
            total_time=2000,
            output_type='timeseries')

        plot_single(
            output,
            variables=plot_variables_list_stochastic,
            out_dir=os.path.join(out_dir, 'stochastic'),
            filename='variables')

    if args.deterministic_divide:
        output = simulate_bioscrape_cobra(
            division=True,
            total_time=6000,
            output_type='unitless')

        plot_multigen(
            output,
            out_dir=os.path.join(out_dir, 'deterministic_divide'),
            filename='division_multigen')

    if args.stochastic_divide:
        output = simulate_bioscrape_cobra(
            stochastic=True,
            division=True,
            total_time=6000,
            external_volume=1e-9*units.L,
            divide_threshold=1500*units.fg,
            output_type='unitless')

        plot_multigen(
            output,
            out_dir=os.path.join(out_dir, 'stochastic_divide'),
            filename='division_multigen')

    if args.deterministic_spatial:


        output = simulate_bioscrape_cobra(
            division=True,
            spatial=True,
            total_time=6000,
            output_type='unitless')

        deterministic_spatial_out_dir = os.path.join(out_dir, 'deterministic_spatial')
        plot_multigen(
            output,
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_multigen',
        )

        plot_fields_snapshots(
            output,
            bounds=BOUNDS,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_snapshots',
        )

        plot_fields_tags(
                        output,
            bounds=BOUNDS,
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            out_dir=deterministic_spatial_out_dir,
            filename='spatial_tags')

    if args.stochastic_spatial:

        output = simulate_bioscrape_cobra(
            stochastic=True,
            division=True,
            spatial=True,
            total_time=6000,
            output_type='unitless')

        stochastic_spatial_out_dir = os.path.join(out_dir, 'stochastic_spatial')
        plot_multigen(
            output,
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_multigen',
        )

        plot_fields_snapshots(
            output,
            bounds=BOUNDS,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_snapshots',
        )

        plot_fields_tags(
            output,
            bounds=BOUNDS,
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            out_dir=stochastic_spatial_out_dir,
            filename='spatial_tags')


if __name__ == '__main__':
    main()
