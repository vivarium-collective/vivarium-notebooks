import os
import sys
import numpy as np
import argparse

# vivarium imports
from vivarium import TreeMass, DivideCondition, MetaDivision
from vivarium.core.process import Composite
from vivarium.library.units import units
from vivarium.core.composition import COMPOSITE_OUT_DIR, FACTORY_KEY, compose_experiment, compartment_in_experiment
from vivarium.core.experiment import Experiment

# vivarium-core processes
from vivarium.processes.clock import Clock

# vivarium-bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# vivarium-cobra imports
from vivarium_cobra import Volume
from vivarium_cobra.processes.configurations import get_iAF1260b_config
from vivarium_cobra.processes.dynamic_fba import DynamicFBA
from vivarium_cobra.processes.local_field import LocalField

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import Lattice, make_lattice_config

# local import
from biocobra.processes.flux_adaptor import FluxAdaptor, DilutionFluxAdaptor, AverageFluxAdaptor
from biocobra.processes.biomass_adaptor import mass_to_concentration, mass_to_count

# plots
from vivarium.plots.simulation_output import plot_simulation_output, plot_variables
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data,
    plot_snapshots,
)
from vivarium_multibody.plots.snapshots import plot_tags


NAME = 'BioscrapeCOBRA'
SBML_FILE_DETERMINISTIC = 'lac_operon/LacOperon_deterministic.xml'
SBML_FILE_STOCHASTIC = 'lac_operon/LacOperon_stochastic.xml'


#choose the SBML file and set other bioscrape parameters
deterministic_bioscrape_config = {
            'sbml_file': SBML_FILE_DETERMINISTIC,
            'stochastic': False,
            'initial_volume': 1,
            'internal_dt': 0.01,
}

stochastic_bioscrape_config = {
            'sbml_file': SBML_FILE_STOCHASTIC,
            'stochastic': True,
            'safe_mode': False,
            'initial_volume': 1,
            'internal_dt': 0.1,
}

# set cobra constrained reactions config
cobra_config = get_iAF1260b_config()

#set up the config for the FluxAdaptor
flux_config = {
    'flux_keys': {
        'Lactose_consumed': {}, #No options specified
        'Glucose_internal': {},  #No options specified
    },
}

dilution_rate_flux_config = {
    'flux_keys': {
        'biomass': {
            'input_type': 'amount'
        }
    }
}

#Here we override the default ports schema of the Biomass species and the k_dilution rate in Bioscrape.
#This is done so they can be set by the Derivers connected to mass and mass flux from Cobra.
schema_override = {
    'bioscrape': {
        'species': {
            'Biomass': {
                '_updater': 'set'  #override bioscrape ('species', 'Biomass') with a 'set' updater
            }
        },
        'rates':{
            'k_dilution__': {
                '_emit': True,  #k_dilution should be emitted so it can be plotted
                '_updater':'set',
            }
        }
    }
}


class BioscrapeCOBRA(Composite):
    defaults = {
        'bioscrape_deterministic': deterministic_bioscrape_config,
        'bioscrape_stochastic': stochastic_bioscrape_config,
        'cobra': cobra_config,
        'flux_adaptor': flux_config,
        'dilution_rate_flux': dilution_rate_flux_config,
        'divide_on': False,  # is division turned on?
        'agent_id': np.random.randint(0, 100),
        'divide_condition': {
            'threshold': 2000 * units.fg},
        'boundary_path': ('boundary',),
        'agents_path': ('agents',),
        'fields_path': ('fields',),
        'dimensions_path': ('dimensions',),
        'daughter_path': tuple(),
        '_schema': schema_override,
        'stochastic': False,  # Is the CRN stochastic or deterministic?
        'spatial_on': False,  # are spatial dynamics used?
        'bioscrape_timestep': 10,
        'cobra_timestep': 10,
        'clock': {
            'time_step': 1.0}
    }

    def __init__(self, config=None):
        super().__init__(config)
        self.config['dilution_rate_flux']['time_step'] = self.config['cobra_timestep']
        self.config['flux_adaptor']['time_step'] = self.config['bioscrape_timestep']

    def generate_processes(self, config):
        processes = {
            'cobra': DynamicFBA(config['cobra']),
            'mass_deriver': TreeMass(),
            'volume_deriver': Volume(),
            'clock': Clock(config['clock']),
        }

        # Process Logic for different kinds of simulations

        # Deterministic case
        if not config['stochastic']:
            # create a deterministic bioscrape model
            processes['bioscrape'] = Bioscrape(config['bioscrape_deterministic'])

            # deterministic simulations have a variable dilution rate
            processes['dilution_rate_adaptor'] = DilutionFluxAdaptor(config["dilution_rate_flux"])

            # flux is computed as an instaneous flux
            processes['flux_adaptor'] = FluxAdaptor(config['flux_adaptor'])

            # biomass is converted to a concentration
            processes['biomass_adaptor'] = mass_to_concentration()

        # Stochastic Case
        else:
            # create a stochastic bioscrape model
            processes['bioscrape'] = Bioscrape(config['bioscrape_stochastic'])

            # flux is computed as an average flux
            processes['flux_adaptor'] = AverageFluxAdaptor(config['flux_adaptor'])

            # biomass is converted to a molecular count
            processes['biomass_adaptor'] = mass_to_count()

        # Division Logic
        if config['divide_on']:
            # division config
            daughter_path = config['daughter_path']
            agent_id = config['agent_id']
            division_config = dict(
                config.get('division', {}),
                daughter_path=daughter_path,
                agent_id=agent_id,
                generator=self)

            processes.update({
                'divide_condition': DivideCondition(config['divide_condition']),
                'division': MetaDivision(division_config)
            })

        # Spatial logic
        if config["spatial_on"]:
            processes.update({'local_field': LocalField()})

        return processes

    def generate_topology(self, config):
        agents_path = config['agents_path']
        fields_path = config['fields_path']
        dimensions_path = config['dimensions_path']
        boundary_path = config['boundary_path']

        topology = {
            'bioscrape': {
                # all species go to a species store on the base level,
                # except Biomass, which goes to the 'globals' store, with variable 'biomass'
                'species': {
                    '_path': ('species',),
                    'Biomass': ('..',) + boundary_path + ('biomass',),
                },
                'delta_species': ('delta_species',),
                'rates': {
                    '_path': ('rates',),
                },
                'globals': boundary_path,
            },
            'cobra': {
                'internal_counts': ('internal_counts',),
                'external': boundary_path + ('external',),
                'exchanges': boundary_path + ('exchange',),
                'reactions': ('reactions',),
                'flux_bounds': ('flux_bounds',),
                'global': boundary_path,
            },
            'flux_adaptor': {
                'inputs': ('delta_species',),
                # 'amounts': boundary_path,
                # connect Bioscrape deltas 'Lactose_consumed' and 'Glucose_internal'
                # to COBRA flux bounds 'EX_lac__D_e' and 'EX_glc__D_e'

                'fluxes': {
                    '_path': ('flux_bounds',),
                    'Lactose_consumed': ('EX_lac__D_e',),
                    'Glucose_internal': ('EX_glc__D_e',),
                }
            },

            'mass_deriver': {
                'global': boundary_path,
            },
            'volume_deriver': {
                'global': boundary_path,
            },
            'biomass_adaptor': {
                'input': boundary_path,
                'output': boundary_path,
            },
            'clock': {
                'global_time': boundary_path + ('time',)
            }
        }

        # Ports added only in the deterministic case
        if not config['stochastic']:
            # Create port biomass flux to the dilution rate computed by the dilution_rate_adaptor process
            topology['dilution_rate_adaptor'] = {
                'inputs': boundary_path,
                'fluxes': {
                    '_path': ('rates',),
                    'biomass': ('k_dilution__',)
                }
            }

        # Ports added only in the stochastic case
        else:
            pass

        if config['divide_on']:
            # connect divide_condition to the mass variable
            topology.update({
                'divide_condition': {
                    'variable': boundary_path + ('mass',),
                    'divide': boundary_path + ('divide',),
                },
                'division': {
                    'global': boundary_path,
                    'agents': agents_path,
                },
            })

        # Ports to use in the spatial case
        if config["spatial_on"]:
            topology.update({'local_field': {
                'exchanges': boundary_path + ('exchange',),
                'location': boundary_path + ('location',),
                'fields': fields_path,
                'dimensions': dimensions_path,
            }})

        return topology

# tests

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
BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 10
#Turn on spatial dynamics
spatial_config = dict(divide_config)
spatial_config["spatial_on"] = True



def test_bioscrape_cobra(total_time=1000):

    bioscrape_composer = BioscrapeCOBRA({})

    initial_state = bioscrape_composer.initial_state()
    initial_state['species']['Glucose_external'] = 10 ** 6
    initial_state['species']['Lactose_external'] = 10 ** 6

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

    bioscrape_composer = BioscrapeCOBRA({'stochastic':True})

    initial_state = bioscrape_composer.initial_state()
    initial_state['species']['Glucose_external'] = 10 ** 6
    initial_state['species']['Lactose_external'] = 10 ** 6

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


def test_bioscrape_cobra_divide():
    total_time = 5000

    division_composite = BioscrapeCOBRA(divide_config)

    # initial state
    initial_state = division_composite.initial_state()
    initial_state['species']['Glucose_external'] = 1e6
    initial_state['species']['Lactose_external'] = 1e6
    initial_state = {
        'agents': {
            agent_id: initial_state}}

    # run simulation
    settings = {
        'outer_path': outer_path,
        'initial_state': initial_state,
        'experiment_id': 'division'}
    division_experiment = compartment_in_experiment(
        division_composite,
        settings=settings,
        initial_state=initial_state)

    # run the experiment and extract the data
    division_experiment.update(total_time)
    division_output = division_experiment.emitter.get_data_unitless()

    final_agents = division_output[total_time]['agents'].keys()
    assert len(final_agents) > 1, 'bioscrapeCOBRA agent did not successfully divide'
    return division_output


def test_bioscrape_cobra_lattice(total_time=100):

    # get initial state
    fields_composer = BioscrapeCOBRA(spatial_config)
    initial_state = fields_composer.initial_state()

    # initial external
    initial_field_concs = {}  #initial_state['boundary']['external']
    initial_field_concs.update({
        'glc__D_e': 10,
        'lcts_e': 10
    })

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
        FACTORY_KEY: {
            'type': Lattice,
            'config': lattice_config},
        'agents': {
            agent_id: {
                FACTORY_KEY: {
                    'type': BioscrapeCOBRA,
                    'config': divide_config
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
    out_dir = os.path.join(COMPOSITE_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='bioscrape_cobra')
    parser.add_argument('--biocobra', '-b', action='store_true', default=False)
    parser.add_argument('--divide', '-d', action='store_true', default=False)
    parser.add_argument('--stochastic', '-s', action='store_true', default=False)
    parser.add_argument('--environment', '-e', action='store_true', default=False)
    parser.add_argument('--all', '-a', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    if args.environment or args.all:
        output = test_bioscrape_cobra_lattice()

        # multigen plots
        plot_settings = {
            'skip_paths': [
                ('external',),
                ('internal_counts',),
            ],
            'remove_zeros': True}
        plot_agents_multigen(
            output, plot_settings, out_dir, 'spatial_multigen')

        agents, fields = format_snapshot_data(output)
        plot_snapshots(
            bounds=BOUNDS,
            agents=agents,
            fields=fields,
            include_fields=['glc__D_e', 'lcts_e'],
            out_dir=out_dir,
            filename='spatial_snapshots')

        tags_data = {'agents': agents, 'fields': fields, 'config': {'bounds': BOUNDS}}
        tags_config = {
            'tagged_molecules': [
                ('species', 'protein_Lactose_Permease',),
            ],
            'out_dir': out_dir,
            'filename': 'spatial_tags'}
        plot_tags(
            data=tags_data,
            plot_config=tags_config
        )


    if args.divide or args.all:
        output = test_bioscrape_cobra_divide()
        # multigen plots
        plot_settings = {
            'skip_paths': [
                ('external',),
                ('internal_counts',),
            ],
            'remove_zeros': True}
        plot_agents_multigen(
            output, plot_settings, out_dir, 'division_multigen')


    if args.stochastic or args.all:
        output = test_bioscrape_cobra_stochastic()

        # plot output
        variables_plot_config = {
            'filename': 'composite_alone_variables',
            'row_height': 2,
            'row_padding': 0.2,
            'column_width': 10,
            'out_dir': out_dir,
            'variables': [
                ('species', 'Glucose_external'),
                ('species', 'Lactose_external'),
                ('species', 'rna_M'),
                ('species', 'protein_betaGal'),
                ('species', 'protein_Lactose_Permease')]}

        plot_variables(output, **variables_plot_config)
        plot_simulation_output(output,
                               out_dir=out_dir,
                               filename='composite_alone_stochastic',
                               )

    if args.biocobra or args.all:
        output = test_bioscrape_cobra()

        # plot output
        variables_plot_config = {
            'filename': 'composite_alone_variables',
            'row_height': 2,
            'row_padding': 0.2,
            'column_width': 10,
            'out_dir': out_dir,
            'variables': [
                ('species', 'Glucose_external'),
                ('species', 'Lactose_external'),
                ('species', 'rna_M'),
                ('species', 'protein_betaGal'),
                ('species', 'protein_Lactose_Permease')]}

        plot_variables(output, **variables_plot_config)
        plot_simulation_output(output,
                               out_dir=out_dir,
                               filename='composite_alone',
                               )

if __name__ == '__main__':
    run_bioscrape_cobra()
