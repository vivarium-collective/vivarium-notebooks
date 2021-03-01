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

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import (
    Lattice, make_lattice_config)

# local import
from bioscrape_cobra.bioscrape_cobra_stochastic import (
    BioscrapeCOBRAstochastic, GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL)
from bioscrape_cobra.bioscrape_cobra_deterministic import BioscrapeCOBRAdeterministic

# plotting
from bioscrape_cobra.plot import (
    plot_multigen, plot_single, plot_fields)

# default variables, which can be varied by simulate_bioscrape_cobra
DEFAULT_EXTERNAL_VOLUME = 1e-13 * units.L
DEFAULT_DIVIDE_THRESHOLD = 2000 * units.fg
INITIAL_GLC = 1e0  # mmolar
INITIAL_LAC = 1e0  # mmolar
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 20

# fixed global variables
COBRA_TIMESTEP = 10
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
                'bin_volume': external_volume}}

    return config



def simulate_bioscrape_cobra(
        division=False,
        stochastic=False,
        spatial=False,
        initial_glucose=1e0,
        initial_lactose=1e0,
        initial_state=None,
        bounds=[20, 20],
        n_bins=[10, 10],
        depth=20,
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
        state = biocobra_composite.initial_state()
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

        plot_fields(
            output,
            bounds=BOUNDS,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            out_dir=deterministic_spatial_out_dir,
            filename='spatial',
        )

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

        plot_fields(
            output,
            bounds=BOUNDS,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            tagged_molecules=[('species', 'protein_Lactose_Permease',)],
            out_dir=stochastic_spatial_out_dir,
            filename='spatial',
        )


if __name__ == '__main__':
    main()
