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
from vivarium.plots.simulation_output import (
    plot_simulation_output, plot_variables)
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data, plot_snapshots)
from vivarium_multibody.plots.snapshots import plot_tags

# default variables, which can be varied by simulate_bioscrape_cobra
DEFAULT_EXTERNAL_VOLUME = 1e-12 * units.L
DEFAULT_DIVIDE_THRESHOLD = 2000 * units.fg

# global variables
GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
SBML_FILE_DETERMINISTIC = 'bioscrape_cobra/LacOperon_deterministic.xml'
SBML_FILE_STOCHASTIC = 'bioscrape_cobra/LacOperon_stochastic.xml'
COBRA_TIMESTEP = 10
BIOSCRAPE_TIMESTEP = 10

# divide config
AGENT_ID = '1'
divide_config = {
    'divide_on': True,
    'AGENT_ID': AGENT_ID,
    'agents_path': ('..', '..', 'agents',),
    'fields_path': ('..', '..', 'fields',),
    'dimensions_path': ('..', '..', 'dimensions',),
    'local_fields': {}}

# spatial config
spatial_config = dict(divide_config)
spatial_config['fields_on'] = True

# lattice config
INITIAL_GLC = 1e0
INITIAL_LAC = 1e0
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 20

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
    ('boundary', 'external', LACTOSE_EXTERNAL)]
plot_variables_list_deterministic.extend(plot_variables_list)

plot_variables_list_stochastic = [
    ('species', GLUCOSE_EXTERNAL),
    ('species', LACTOSE_EXTERNAL)]
plot_variables_list_stochastic.extend(plot_variables_list)


# helper functions
def get_bioscrape_cobra_config(
        division=False,
        divide_threshold=DEFAULT_DIVIDE_THRESHOLD,
        external_volume=DEFAULT_EXTERNAL_VOLUME,
):
    """ create a generic config dict for bioscrape_cobra composers """
    config = {
        'local_fields': {
            'bin_volume': external_volume}}
    if division:
        config['divide_condition'] = {
            'threshold': divide_threshold}
    return config


def put_bioscrape_cobra_in_lattice(
        biocobra_composer,
        biocobra_config,
        field_concentrations,
        bounds=BOUNDS,
        n_bins=NBINS,
        depth=DEPTH,
        agent_ids=[AGENT_ID],
):
    """ configure lattice compartment
    :return: a hierarchy dict for compose_experiment to initialize
    """
    lattice_config_kwargs = {
        'bounds': bounds,
        'n_bins': n_bins,
        'depth': depth,
        'concentrations': field_concentrations,
        'diffusion': 1e-1,
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
                    'config': biocobra_config}}
            for agent_id in agent_ids}}

    return hierarchy


def simulate_bioscrape_cobra(
        division=False,
        stochastic=False,
        spatial=False,
        initial_glucose=1e0,
        initial_lactose=1e0,
        divide_threshold=2000*units.fg,
        agent_ids=None,
        total_time=100,
        output_type=None,
):
    """ Simulation function for BioscrapeCOBRA """

    # get the composer and configuration
    if agent_ids is None:
        agent_ids = [AGENT_ID]
    if stochastic:
        biocobra_composer = BioscrapeCOBRAstochastic
        biocobra_config = get_bioscrape_cobra_config(
            division=division,
            divide_threshold=divide_threshold)

        # get initial state from composer
        composer_instance = biocobra_composer(biocobra_config)
        initial_state = composer_instance.initial_state()

    else:
        biocobra_composer = BioscrapeCOBRAdeterministic
        biocobra_config = get_bioscrape_cobra_config(
            division=division,
            divide_threshold=divide_threshold)

        # get initial state from composer
        composer_instance = biocobra_composer(biocobra_config)
        initial_state = composer_instance.initial_state()
        initial_state['boundary']['external'] = {
            GLUCOSE_EXTERNAL: initial_glucose,
            LACTOSE_EXTERNAL: initial_lactose}

    if spatial:
        # place the biocobra_composer in a hierarchy with lattice
        hierarchy = put_bioscrape_cobra_in_lattice(
            biocobra_composer, biocobra_config)

    elif division:
        # place the biocobra_composer in an embedded hierarchy under ('agents', AGENT_ID)
        hierarchy = {
            'agents': {
                agent_id: {
                    COMPOSER_KEY: {
                        'type': biocobra_composer,
                        'config': biocobra_config}}
                for agent_id in agent_ids}}
    else:
        hierarchy = {
            COMPOSER_KEY: {
                'type': biocobra_composer,
                'config': biocobra_config}}

    # make experiment with helper function compose_experiment()
    experiment_settings = {
        'initial_state': initial_state,
        'experiment_id': f"{division} {stochastic} {spatial}"}
    biocobra_experiment = compose_experiment(
        hierarchy=hierarchy,
        settings=experiment_settings)

    # run the experiment
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
    parser.add_argument('--deterministic', '-d', action='store_true', default=False)
    parser.add_argument('--stochastic', '-s', action='store_true', default=False)
    args = parser.parse_args()

    if args.deterministic:
        biocobra_out_dir = os.path.join(out_dir, 'deterministic')
        output = simulate_bioscrape_cobra(
            total_time=100,
            output_type='timeseries')

        # plot output
        variables_plot_config = {
            'out_dir': biocobra_out_dir, 'filename': 'variables',
            'row_height': 2, 'row_padding': 0.2, 'column_width': 10,
            'variables': plot_variables_list_deterministic}
        plot_variables(
            output, **variables_plot_config)

    if args.stochastic:
        pass


if __name__ == '__main__':
    main()
