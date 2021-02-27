"""
==============================================
Simulation helper functions for BioscrapeCOBRA
==============================================
"""
import os
import argparse

# vivarium imports
from vivarium.core.composition import compose_experiment, COMPOSER_KEY
from vivarium.library.units import units

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import (
    Lattice, make_lattice_config)

# local import
from bioscrape_cobra.bioscrape_cobra_stochastic import BioscrapeCOBRAstochastic, GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL
from bioscrape_cobra.bioscrape_cobra_deterministic import BioscrapeCOBRAdeterministic


# global variables
GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
SBML_FILE_DETERMINISTIC = 'bioscrape_cobra/LacOperon_deterministic.xml'
COBRA_TIMESTEP = 10
BIOSCRAPE_TIMESTEP = 10

# divide config
agent_id = '1'
outer_path = ('agents', agent_id,)
divide_config = {
    'divide_on': True,
    'agent_id': agent_id,
    'agents_path': ('..', '..', 'agents',),
    'fields_path': ('..', '..', 'fields',),
    'dimensions_path': ('..', '..', 'dimensions',),
    'local_fields': {}}

# spatial config
spatial_config = dict(divide_config)
spatial_config['fields_on'] = True

# lattice environment spatial config
INITIAL_GLC = 1e0
INITIAL_LAC = 1e0
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 20


def get_deterministic_composer(
        division=False,
        divide_threshold=2000*units.fg,
        external_volume=1e-12 * units.L,
):
    biocobra_config = {
        'local_fields': {
            'bin_volume': external_volume}}
    if division:
        biocobra_config['divide_condition'] = {'threshold': divide_threshold}

    return BioscrapeCOBRAdeterministic(biocobra_config)


def get_stochastic_composer(
        division=False,
        external_volume=1e-12 * units.L,
):
    biocobra_config = {
        'local_fields': {
            'bin_volume': external_volume}}
    return BioscrapeCOBRAstochastic(biocobra_config)


def put_in_lattice(
        field_concentrations,
        lattice_timestep=COBRA_TIMESTEP,
        bounds=BOUNDS,
        n_bins=NBINS,
        depth=DEPTH,
):
    # configure lattice compartment
    lattice_config_kwargs = {
        'bounds': bounds,
        'n_bins': n_bins,
        'depth': depth,
        'concentrations': field_concentrations,
        'diffusion': 1e-1,
        'time_step': lattice_timestep}
    lattice_config = make_lattice_config(**lattice_config_kwargs)

    # declare the hierarchy
    hierarchy = {
        COMPOSER_KEY: {
            'type': Lattice,
            'config': lattice_config},
        'agents': {
            agent_id: {
                COMPOSER_KEY: {
                    'type': BioscrapeCOBRAdeterministic,
                    'config': spatial_config}}}}

    return hierarchy



def simulate_bioscrape_cobra(
        division=False,
        stochastic=False,
        spatial=False,
        initial_glucose=1e0,
        initial_lactose=1e0,
        divide_threshold=2000*units.fg,
        total_time=100,
):
    # get the composite composer
    if stochastic:
        composer = get_stochastic_composer()
        initial_state = composer.initial_state()
    else:
        composer = get_deterministic_composer()

        # get initial state
        initial_state = composer.initial_state()
        # initial_state['boundary']['biomass'] = 0.00182659297 * units.mmolar  # TODO -- is this still needed?
        initial_state['boundary']['external'] = {
            GLUCOSE_EXTERNAL: initial_glucose,
            LACTOSE_EXTERNAL: initial_lactose}


    if spatial:
        # place it in a hierarchy
        hierarchy = put_in_lattice(composer)
    else:
        hierarchy = composite

    # make experiment with helper function compose_experiment()
    experiment_settings = {
        'initial_state': initial_state,
        'experiment_id': f"{division} {stochastic} {spatial}"}
    spatial_experiment = compose_experiment(
        hierarchy=hierarchy,
        settings=experiment_settings)

    # run the experiment
    spatial_experiment.update(total_time)

    # retrieve the data
    data = spatial_experiment.emitter.get_data_unitless()

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
    args = parser.parse_args()

    if args.deterministic:
        biocobra_out_dir = os.path.join(out_dir, 'single')
        simulate_bioscrape_cobra(
            total_time=100
        )


if __name__ == '__main__':
    main()
