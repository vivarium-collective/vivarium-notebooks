"""
==============================================
Simulation helper functions for BioscrapeCOBRA
==============================================
"""
import os
import argparse
import copy

# vivarium imports
from vivarium.core.composition import compose_experiment, COMPOSER_KEY
from vivarium.library.units import units

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import (
    Lattice, make_lattice_config)

# local import
from bioscrape_cobra.bioscrape_cobra_stochastic import BioscrapeCOBRAstochastic, GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL
from bioscrape_cobra.bioscrape_cobra_deterministic import BioscrapeCOBRAdeterministic

# plotting
from vivarium.plots.simulation_output import plot_variables
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data, plot_snapshots)
from vivarium_multibody.plots.snapshots import plot_tags

# default variables, which can be varied by simulate_bioscrape_cobra
DEFAULT_EXTERNAL_VOLUME = 1e-13 * units.L
DEFAULT_DIVIDE_THRESHOLD = 2000 * units.fg
INITIAL_GLC = 1e0
INITIAL_LAC = 1e0
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 20

# fixed global variables
GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
SBML_FILE_DETERMINISTIC = 'bioscrape_cobra/LacOperon_deterministic.xml'
SBML_FILE_STOCHASTIC = 'bioscrape_cobra/LacOperon_stochastic.xml'
COBRA_TIMESTEP = 10
BIOSCRAPE_TIMESTEP = 10

# divide config
INITIAL_AGENT_ID = '1'
divide_config = {
    'divide_on': True,
    'agent_id': INITIAL_AGENT_ID,
    'agents_path': ('..', '..', 'agents',),
    'fields_path': ('..', '..', 'fields',),
    'dimensions_path': ('..', '..', 'dimensions',),
    'local_fields': {}}

# spatial config
spatial_config = dict(divide_config)
spatial_config['fields_on'] = True

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


# helper functions
def get_bioscrape_cobra_config(
        spatial=False,
        division=False,
        divide_threshold=DEFAULT_DIVIDE_THRESHOLD,
):
    """ create a generic config dict for bioscrape_cobra composers """
    agent_id = INITIAL_AGENT_ID
    external_volume = DEFAULT_EXTERNAL_VOLUME

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


def put_bioscrape_cobra_in_lattice(
        biocobra_composer,
        biocobra_config,
        field_concentrations=None,
        diffusion_rate=1e-1,
        bounds=BOUNDS,
        n_bins=NBINS,
        depth=DEPTH,
        agent_ids=None,
):
    """ configure lattice compartment
    :return: a hierarchy dict for compose_experiment to initialize
    """
    if field_concentrations is None:
        field_concentrations = {
            GLUCOSE_EXTERNAL: INITIAL_GLC,
            LACTOSE_EXTERNAL: INITIAL_LAC}
    agent_id = INITIAL_AGENT_ID
    lattice_config_kwargs = {
        'bounds': bounds,
        'n_bins': n_bins,
        'depth': depth,
        'concentrations': field_concentrations,
        'diffusion': diffusion_rate,
        'time_step': COBRA_TIMESTEP}
    lattice_config = make_lattice_config(**lattice_config_kwargs)

    # declare the hierarchy
    hierarchy = {
        COMPOSER_KEY: {
            'type': Lattice,
            'config': lattice_config},
        'agents': {
            agent_id: {
                COMPOSER_KEY: {
                    'type': biocobra_composer,
                    'config': biocobra_config}}}}
    return hierarchy


def simulate_bioscrape_cobra(
        division=False,
        stochastic=False,
        spatial=False,
        initial_glucose=1e0,
        initial_lactose=1e0,
        diffusion_rate=1e-1,
        divide_threshold=2000*units.fg,
        halt_threshold=32,
        total_time=100,
        output_type=None,
):
    """ Simulation function for BioscrapeCOBRA """
    agent_id = INITIAL_AGENT_ID

    # get the BioscrapeCOBRA composer
    if stochastic:
        biocobra_composer = BioscrapeCOBRAstochastic
    else:
        biocobra_composer = BioscrapeCOBRAdeterministic

    # make the config
    biocobra_config = get_bioscrape_cobra_config(
        spatial=spatial,
        division=division,
        divide_threshold=divide_threshold)

    # initial state from composer
    composer_instance = biocobra_composer(biocobra_config)
    initial_state = composer_instance.initial_state()
    initial_state['boundary']['external'] = {
        GLUCOSE_EXTERNAL: initial_glucose,
        LACTOSE_EXTERNAL: initial_lactose}

    # make the hierarchy
    if spatial:
        field_concentrations = {
            GLUCOSE_EXTERNAL: initial_glucose,
            LACTOSE_EXTERNAL: initial_lactose}

        # spatial places the agent in a hierarchy with a Lattice composite in the environment
        hierarchy = put_bioscrape_cobra_in_lattice(
            biocobra_composer=biocobra_composer,
            biocobra_config=biocobra_config,
            diffusion_rate=diffusion_rate,
            field_concentrations=field_concentrations)

        initial_state = {
            'agents': {
                agent_id: initial_state}}

    elif division:
        # division requires the agent to be embedded in a hierarchy
        hierarchy = {
            'agents': {
                agent_id: {
                    COMPOSER_KEY: {
                        'type': biocobra_composer,
                        'config': biocobra_config}}}}
        initial_state = {
            'agents': {
                agent_id: initial_state}}

    else:
        hierarchy = {
            COMPOSER_KEY: {
                'type': biocobra_composer,
                'config': biocobra_config}}

    # make the experiment with helper function compose_experiment
    experiment_settings = {
        'initial_state': initial_state,
        # 'display_info': False,
        'experiment_id': f"{'stochastic' if stochastic else 'deterministic'}_"
                         f"{'division' if division else ''}_"
                         f"{'spatial' if spatial else ''}"}
    biocobra_experiment = compose_experiment(
        hierarchy=hierarchy,
        settings=experiment_settings)

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



def plot_bioscrape_cobra(
    stochastic=False,
):
    # plotting
    if stochastic:
        pass
    else:
        pass


def main():
    out_dir = os.path.join(
        'out', 'bioscrape_cobra')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='bioscrape_cobra')
    parser.add_argument('--deterministic', '-a', action='store_true', default=False)
    parser.add_argument('--stochastic', '-b', action='store_true', default=False)
    parser.add_argument('--deterministic_divide', '-c', action='store_true', default=False)
    parser.add_argument('--stochastic_divide', '-d', action='store_true', default=False)
    parser.add_argument('--deterministic_spatial', '-e', action='store_true', default=False)
    parser.add_argument('--stochastic_spatial', '-f', action='store_true', default=False)
    args = parser.parse_args()

    if args.deterministic:
        biocobra_out_dir = os.path.join(out_dir, 'deterministic')
        output = simulate_bioscrape_cobra(
            total_time=2000,
            output_type='timeseries')

        # plot output
        variables_plot_config = {
            'out_dir': biocobra_out_dir, 'filename': 'variables',
            'row_height': 2, 'row_padding': 0.2, 'column_width': 10,
            'variables': plot_variables_list_deterministic}
        plot_variables(
            output, **variables_plot_config)

    if args.stochastic:
        biocobra_out_dir = os.path.join(out_dir, 'stochastic')
        output = simulate_bioscrape_cobra(
            stochastic=True,
            total_time=2000,
            output_type='timeseries')

        # plot output
        variables_plot_config = {
            'out_dir': biocobra_out_dir, 'filename': 'variables',
            'row_height': 2, 'row_padding': 0.2, 'column_width': 10,
            'variables': plot_variables_list_stochastic}
        plot_variables(
            output, **variables_plot_config)

    if args.deterministic_divide:
        biocobra_out_dir = os.path.join(out_dir, 'deterministic_divide')
        output = simulate_bioscrape_cobra(
            division=True,
            total_time=6000,
            output_type='unitless')

        # multigen plot
        plot_settings = {
            'skip_paths': [
                ('internal_counts',),
                ('cobra_external',)],
            'remove_zeros': True}
        plot_agents_multigen(
            output,
            plot_settings,
            biocobra_out_dir,
            'division_multigen')

    if args.stochastic_divide:
        biocobra_out_dir = os.path.join(out_dir, 'stochastic_divide')
        output = simulate_bioscrape_cobra(
            stochastic=True,
            division=True,
            total_time=6000,
            output_type='unitless')

        # multigen plot
        plot_settings = {
            'skip_paths': [
                ('internal_counts',),
                ('cobra_external',)],
            'remove_zeros': False}
        plot_agents_multigen(
            output,
            plot_settings,
            biocobra_out_dir,
            'division_multigen')

    if args.deterministic_spatial:

        biocobra_out_dir = os.path.join(out_dir, 'deterministic_spatial')
        output = simulate_bioscrape_cobra(
            division=True,
            spatial=True,
            total_time=6000,
            output_type='unitless')

        # multigen plots
        plot_settings = {
            'skip_paths': [
                ('internal_counts',),
                ('cobra_external',)],
            'n_snapshots': 5,
            'remove_zeros': True}
        plot_agents_multigen(
            output, plot_settings, biocobra_out_dir, 'spatial_multigen')

        agents, fields = format_snapshot_data(output)
        plot_snapshots(
            bounds=BOUNDS,
            agents=agents,
            fields=fields,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            out_dir=biocobra_out_dir,
            filename='spatial_snapshots')

        tags_data = {
            'agents': agents,
            'fields': fields,
            'config': {'bounds': BOUNDS}}
        tags_config = {
            'tagged_molecules': [
                ('species', 'protein_Lactose_Permease',)],
            'n_snapshots': 5,
            'out_dir': biocobra_out_dir,
            'filename': 'spatial_tags'}
        plot_tags(
            data=tags_data,
            plot_config=tags_config)

    if args.stochastic_spatial:
        biocobra_out_dir = os.path.join(out_dir, 'stochastic_spatial')
        output = simulate_bioscrape_cobra(
            stochastic=True,
            division=True,
            spatial=True,
            total_time=10000,
            output_type='unitless')

        # multigen plots
        plot_settings = {
            'skip_paths': [
                ('internal_counts',),
                ('cobra_external',)],
            'n_snapshots': 5,
            'remove_zeros': True}
        plot_agents_multigen(
            output, plot_settings, biocobra_out_dir, 'spatial_multigen')

        agents, fields = format_snapshot_data(output)
        plot_snapshots(
            bounds=BOUNDS,
            agents=agents,
            fields=fields,
            include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
            out_dir=biocobra_out_dir,
            filename='spatial_snapshots')

        tags_data = {
            'agents': agents,
            'fields': fields,
            'config': {'bounds': BOUNDS}}
        tags_config = {
            'tagged_molecules': [
                ('species', 'protein_Lactose_Permease',)],
            'n_snapshots': 5,
            'out_dir': biocobra_out_dir,
            'filename': 'spatial_tags'}
        plot_tags(
            data=tags_data,
            plot_config=tags_config)


if __name__ == '__main__':
    main()
