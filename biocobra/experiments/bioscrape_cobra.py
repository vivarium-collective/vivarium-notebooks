import os
import sys
import argparse

# vivarium imports
from vivarium.core.composition import EXPERIMENT_OUT_DIR, COMPOSER_KEY, compose_experiment, composer_in_experiment
from vivarium.core.experiment import Experiment

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import Lattice, make_lattice_config

# local import
from biocobra.composites.bioscrape_cobra import NAME, GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL, BioscrapeCOBRA
# from biocobra.processes.biomass_adaptor import mass_to_concentration, mass_to_count

# plots
from vivarium.plots.simulation_output import plot_simulation_output, plot_variables
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data,
    plot_snapshots,
)
from vivarium_multibody.plots.snapshots import plot_tags


# division test config
agent_id = '1'
outer_path = ('agents', agent_id,)
divide_config = {
    'divide_on': True,
    'agent_id': agent_id,
    'agents_path': ('..', '..', 'agents',),
    'fields_path': ('..', '..', 'fields',),
    'dimensions_path': ('..', '..', 'dimensions',)}


# lattice environment test config
BOUNDS = [10, 10]
NBINS = [5, 5]
DEPTH = 10
#Turn on spatial dynamics
spatial_config = dict(divide_config)
spatial_config["spatial_on"] = True



def test_bioscrape_cobra(total_time=1000):

    bioscrape_composer = BioscrapeCOBRA({})

    initial_state = bioscrape_composer.initial_state()
    initial_state['boundary']['external'][GLUCOSE_EXTERNAL] = 1e6
    initial_state['boundary']['external'][LACTOSE_EXTERNAL] = 1e5
    # initial_state['species'][GLUCOSE_EXTERNAL] = 1e6
    # initial_state['species'][LACTOSE_EXTERNAL] = 1e5

    # make experiment
    bioscrape_composite = bioscrape_composer.generate()
    bioscrape_experiment = Experiment(
        dict(
            processes=bioscrape_composite['processes'],
            topology=bioscrape_composite['topology'],
            initial_state=initial_state,))

    bioscrape_experiment.update(total_time)
    timeseries = bioscrape_experiment.emitter.get_timeseries()
    return timeseries


def test_bioscrape_cobra_stochastic(total_time=1000):

    stochastic_biocobra_composer = BioscrapeCOBRA({'stochastic': True})

    initial_state = stochastic_biocobra_composer.initial_state()
    initial_state['boundary']['external'][GLUCOSE_EXTERNAL] = 1e6
    initial_state['boundary']['external'][LACTOSE_EXTERNAL] = 1e5
    # initial_state['species'][GLUCOSE_EXTERNAL] = 1e6
    # initial_state['species'][LACTOSE_EXTERNAL] = 1e5

    # make experiment
    stochastic_biocobra_composite = stochastic_biocobra_composer.generate()
    stochastic_biocobra_experiment = Experiment(
        dict(
            processes=stochastic_biocobra_composite['processes'],
            topology=stochastic_biocobra_composite['topology'],
            initial_state=initial_state, ))

    stochastic_biocobra_experiment.update(total_time)
    timeseries = stochastic_biocobra_experiment.emitter.get_timeseries()
    return timeseries


def test_bioscrape_cobra_divide():
    total_time = 2000

    division_composite = BioscrapeCOBRA(divide_config)

    # initial state
    initial_state = division_composite.initial_state()
    initial_state['species'][GLUCOSE_EXTERNAL] = 1e6
    initial_state['species'][LACTOSE_EXTERNAL] = 1e6
    initial_state = {
        'agents': {
            agent_id: initial_state}}

    # run simulation
    settings = {
        'experiment_id': 'division'}
    division_experiment = composer_in_experiment(
        division_composite,
        settings=settings,
        outer_path=outer_path,
        initial_state=initial_state)

    # run the experiment and extract the data
    division_experiment.update(total_time)
    division_output = division_experiment.emitter.get_data_unitless()

    final_agents = division_output[total_time]['agents'].keys()
    assert len(final_agents) > 1, 'bioscrapeCOBRA agent did not successfully divide'
    return division_output


def test_bioscrape_cobra_lattice(total_time=2500):

    # get initial state
    fields_composer = BioscrapeCOBRA(spatial_config)
    initial_state = fields_composer.initial_state()
    # initial_state['species'][GLUCOSE_EXTERNAL] = 1e6  # TODO (Eran): remove this!
    # initial_state['species'][LACTOSE_EXTERNAL] = 1e6  # TODO (Eran): remove this!

    # initial external
    initial_field_concs = {
        GLUCOSE_EXTERNAL: 10,
        LACTOSE_EXTERNAL: 10,
    }
    # initial_field_concs = {}  #initial_state['boundary']['external']
    # initial_field_concs.update({
    #     'glc__D_e': 10,
    #     'lcts_e': 10
    # })

    import ipdb; ipdb.set_trace()
    # TODO -- connect GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL

    # initial agents
    initial_state = {
        'agents': {
            agent_id: initial_state}}

    # configure lattice compartment
    lattice_config_kwargs = {
        'bounds': BOUNDS,
        'n_bins': NBINS,
        'depth': DEPTH,
        'concentrations': initial_field_concs}

    lattice_config = make_lattice_config(**lattice_config_kwargs)

    # declare the hierarchy
    hierarchy = {
        COMPOSER_KEY: {
            'type': Lattice,
            'config': lattice_config},
        'agents': {
            agent_id: {
                COMPOSER_KEY: {
                    'type': BioscrapeCOBRA,
                    'config': spatial_config
                }
            }
        }}

    # make experiment with helper function compose_experiment()
    experiment_settings = {
        'initial_state': initial_state,
        'experiment_id': 'spatial_environment'}
    spatial_experiment = compose_experiment(
        hierarchy=hierarchy,
        settings=experiment_settings)

    spatial_experiment.update(total_time)
    data = spatial_experiment.emitter.get_data_unitless()
    return data


def run_bioscrape_cobra():
    out_dir = EXPERIMENT_OUT_DIR #os.path.join(EXPERIMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='bioscrape_cobra')
    parser.add_argument('--biocobra', '-b', action='store_true', default=False)
    parser.add_argument('--stochastic', '-s', type=int, const=1, nargs='?', default=False)
    parser.add_argument('--divide', '-d', action='store_true', default=False)
    parser.add_argument('--environment', '-e', action='store_true', default=False)
    parser.add_argument('--all', '-a', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    if args.biocobra or args.all:
        biocobra_out_dir = os.path.join(out_dir, 'biocobra')
        output = test_bioscrape_cobra()

        # plot output
        variables_plot_config = {
            'filename': 'composite_alone_variables',
            'row_height': 2,
            'row_padding': 0.2,
            'column_width': 10,
            'out_dir': biocobra_out_dir,
            'variables': [
                ('boundary', 'external', GLUCOSE_EXTERNAL),
                ('boundary', 'external', LACTOSE_EXTERNAL),
                # ('species', GLUCOSE_EXTERNAL),
                # ('species', LACTOSE_EXTERNAL),
                ('species', 'rna_M'),
                ('species', 'protein_betaGal'),
                ('species', 'protein_Lactose_Permease')]}

        plot_variables(output, **variables_plot_config)
        plot_simulation_output(output,
                               out_dir=biocobra_out_dir,
                               filename='composite_alone')

    if args.stochastic or args.all:
        stoch_out_dir = os.path.join(out_dir, 'biocobra_stochastic')
        n_runs = int(args.stochastic) or 1
        for n in range(1, n_runs+1):
            output = test_bioscrape_cobra_stochastic()

            # plot output
            variables_plot_config = {
                'filename': f'composite_alone_variables_{n}',
                'row_height': 2,
                'row_padding': 0.2,
                'column_width': 10,
                'out_dir': stoch_out_dir,
                'variables': [
                    ('boundary', 'external', GLUCOSE_EXTERNAL),
                    ('boundary', 'external', LACTOSE_EXTERNAL),
                    # ('species', GLUCOSE_EXTERNAL),
                    # ('species', LACTOSE_EXTERNAL),
                    ('species', 'rna_M'),
                    ('species', 'protein_betaGal'),
                    ('species', 'protein_Lactose_Permease')]}

            plot_variables(output, **variables_plot_config)
            plot_simulation_output(output,
                                   out_dir=stoch_out_dir,
                                   filename=f'composite_alone_stochastic_{n}')

    if args.divide or args.all:
        div_out_dir = os.path.join(out_dir, 'biocobra_divide')
        output = test_bioscrape_cobra_divide()
        # multigen plots
        plot_settings = {
            'skip_paths': [
                ('external',),
                ('internal_counts',),
            ],
            'remove_zeros': True}
        plot_agents_multigen(
            output, plot_settings, div_out_dir, 'division_multigen')

    if args.environment or args.all:
        env_out_dir = os.path.join(out_dir, 'biocobra_environment')
        output = test_bioscrape_cobra_lattice()

        # multigen plots
        plot_settings = {
            'skip_paths': [
                ('external',),
                ('internal_counts',),
            ],
            'remove_zeros': True}
        plot_agents_multigen(
            output, plot_settings, env_out_dir, 'spatial_multigen')

        agents, fields = format_snapshot_data(output)
        plot_snapshots(
            bounds=BOUNDS,
            agents=agents,
            fields=fields,
            include_fields=['glc__D_e', 'lcts_e'],
            out_dir=env_out_dir,
            filename='spatial_snapshots')

        tags_data = {'agents': agents, 'fields': fields, 'config': {'bounds': BOUNDS}}
        tags_config = {
            'tagged_molecules': [
                ('species', 'protein_Lactose_Permease',),
            ],
            'out_dir': env_out_dir,
            'filename': 'spatial_tags'}
        plot_tags(
            data=tags_data,
            plot_config=tags_config
        )



if __name__ == '__main__':
    run_bioscrape_cobra()
