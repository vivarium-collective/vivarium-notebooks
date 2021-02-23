import os
import numpy as np
import argparse

# vivarium-core processes
from vivarium import (
    TreeMass, Clock, MassToMolar, MassToCount, CountsToMolar,
    DivideCondition, MetaDivision, MolarToCounts, StripUnits)
from vivarium.core.experiment import Experiment
from vivarium.core.process import Composer
from vivarium.library.units import units
from vivarium.core.composition import (
    compose_experiment, EXPERIMENT_OUT_DIR, COMPOSER_KEY)

# vivarium-bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# vivarium-cobra imports
from vivarium_cobra import Volume, LocalField
from vivarium_cobra.processes.configurations import get_iAF1260b_config
from vivarium_cobra.processes.dynamic_fba import DynamicFBA

# vivarium-multibody imports
from vivarium_multibody.composites.lattice import Lattice, make_lattice_config

# local imports
from biocobra.processes.flux_adaptor import AverageFluxAdaptor

# plots
from vivarium.plots.simulation_output import plot_simulation_output, plot_variables
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data,
    plot_snapshots,
)
from vivarium_multibody.plots.snapshots import plot_tags

GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
SBML_FILE_STOCHASTIC = 'lac_operon/LacOperon_stochastic.xml'
COBRA_TIMESTEP = 2
BIOSCRAPE_TIMESTEP = 2

# choose the SBML file and set other bioscrape parameters
stochastic_bioscrape_config = {
    'sbml_file': SBML_FILE_STOCHASTIC,
    'stochastic': True,
    'safe_mode': False,
    'initial_volume': 1,
    'internal_dt': 0.1,
}

# set cobra constrained reactions config
cobra_config = get_iAF1260b_config()

# set up the config for the FluxAdaptor
flux_config = {
    'flux_keys': {
        'Lactose_consumed': {
            'input_type': 'delta',
            'window_size': 5},
        'Glucose_internal': {
            'input_type': 'delta',
            'window_size': 5},
    },
}
mass_mw_config = {
    'molecular_weights': {
        'mass': 1.0 * units.fg / units.molec
    }
}
# convert stochastic delta counts to delta concentrations for constraining FBA
delta_counts_to_concs_config = {
    'keys': [
        'Glucose_internal',
        'Lactose_consumed',
    ]
}
# configures the counts deriver to convert field concentrations to counts for stochastic sims
field_counts_deriver_config = {
    'keys': [GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL]}

# configuration for strip units deriver, which converts and removes specified units
strip_units_config = {
    'keys': [
        'mass', 'volume', 'density',
        'biomass'
    ],
    'convert': {
        # 'biomass': units.mmolar,
        'mass': units.ug,
    }}

# set mass threshold for division
divide_config = {'threshold': 2000 * units.fg}

# Here we override the default ports schema of the Biomass species and the k_dilution rate in Bioscrape.
# This is done so they can be set by the Derivers connected to mass and mass flux from Cobra.
schema_override = {
    'bioscrape': {
        'species': {
            'Biomass': {
                # '_default': 0.00166,
                '_updater': 'set',  # override bioscrape ('species', 'Biomass') with a 'set' updater
            },
            'Glucose_external': {
                '_divider': 'set',
                '_updater': 'null',
            },
            'Lactose_external': {
                '_divider': 'set',
                '_updater': 'null',
            },
            'dna_Lac_Operon': {
                '_divider': 'set',
            },
        },
        'rates': {
            'k_dilution__': {
                '_emit': True,  # k_dilution should be emitted so it can be plotted
                '_updater': 'set',
            }
        }
    }
}


class BioscrapeCOBRAstochastic(Composer):
    defaults = {
        'divide_on': False,  # is division turned on?
        'fields_on': False,  # are spatial dynamics used?
        'bioscrape': stochastic_bioscrape_config,
        'cobra': cobra_config,
        'flux_adaptor': flux_config,
        'mass_to_counts': mass_mw_config,
        'delta_counts_to_concs': delta_counts_to_concs_config,
        'field_counts_deriver': field_counts_deriver_config,
        'strip_units': strip_units_config,
        'local_fields': {},
        'agent_id': np.random.randint(0, 100),
        'divide_condition': divide_config,
        'boundary_path': ('boundary',),
        'agents_path': ('agents',),
        'fields_path': ('fields',),
        'dimensions_path': ('dimensions',),
        'daughter_path': tuple(),
        '_schema': schema_override,
        'bioscrape_timestep': BIOSCRAPE_TIMESTEP,
        'cobra_timestep': COBRA_TIMESTEP,
        'clock': {'time_step': 1.0}}

    def __init__(self, config=None):
        super().__init__(config)
        # configure timesteps

        self.config['bioscrape']['time_step'] = self.config['bioscrape_timestep']
        self.config['flux_adaptor']['time_step'] = self.config['bioscrape_timestep']
        self.config['cobra']['time_step'] = self.config['cobra_timestep']
        self.config['clock']['time_step'] = min(self.config['cobra_timestep'], self.config['bioscrape_timestep'])

        # configure local fields
        if not self.config['fields_on']:
            self.config['local_fields'].update({'nonspatial': True})

    def generate_processes(self, config):
        processes = {
            'cobra': DynamicFBA(config['cobra']),
            'bioscrape': Bioscrape(config['bioscrape']),
            'local_field': LocalField(config['local_fields']),
            'clock': Clock(config['clock']),
            'mass_deriver': TreeMass(),
            'volume_deriver': Volume(),
            'delta_counts_to_concs': CountsToMolar(config['delta_counts_to_concs']),
            'flux_adaptor': AverageFluxAdaptor(config['flux_adaptor']),
            'biomass_adaptor': MassToCount(config['mass_to_counts']),
            'field_counts_deriver': MolarToCounts(config['field_counts_deriver']),
            'strip_units': StripUnits(config['strip_units']),
        }

        # Division Logic
        if config['divide_on']:
            # division config
            daughter_path = config['daughter_path']
            agent_id = config['agent_id']
            division_config = dict(
                config.get('division', {}),
                daughter_path=daughter_path,
                agent_id=agent_id,
                composer=self)

            processes.update({
                'divide_condition': DivideCondition(config['divide_condition']),
                'division': MetaDivision(division_config)})

        return processes

    def generate_topology(self, config):
        agents_path = config['agents_path']
        fields_path = config['fields_path']
        dimensions_path = config['dimensions_path']
        boundary_path = config['boundary_path']
        unitless_boundary_path = boundary_path + ('no_units',)

        topology = {
            'bioscrape': {
                # all species go to a species store on the base level,
                # except Biomass, which goes to the 'boundary' store, with variable 'biomass'
                'species': {
                    '_path': ('species',),
                    'Biomass': ('..',) + unitless_boundary_path + ('biomass',),
                },
                'delta_species': ('delta_species',),
                'rates': ('rates',),
                'globals': unitless_boundary_path,
            },
            'cobra': {
                'internal_counts': ('internal_counts',),
                # 'external': boundary_path + ('external',),
                'external': ('cobra_external',),  # These are handled separately from the external fields
                'exchanges': {
                    # connect only glc__D_e lac__D_e to boundary exchanges that update fields
                    '_path': ('hidden_exchanges',),
                    'glc__D_e': ('..',) + boundary_path + ('exchanges', GLUCOSE_EXTERNAL,),
                    'lac__D_e': ('..',) + boundary_path + ('exchanges', LACTOSE_EXTERNAL,),
                },
                'reactions': ('reactions',),
                'flux_bounds': ('flux_bounds',),
                'global': boundary_path,
            },
            'flux_adaptor': {
                'inputs': ('delta_concentrations',),
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
                'input': {
                    '_path': boundary_path,
                    'mass': ('mass',)
                },
                'output': {
                    '_path': boundary_path,
                    'mass': ('biomass',)
                },
                'global': boundary_path,
            },
            'clock': {
                'global_time': boundary_path + ('time',)
            },
            'strip_units': {
                'units': boundary_path,
                'no_units': unitless_boundary_path,
            },
            'local_field': {
                'exchanges': boundary_path + ('exchanges',),
                'location': boundary_path + ('location',),
                # connect fields directly to external port if fields_on is False
                'fields': fields_path if config['fields_on'] else boundary_path + ('external',),
                'dimensions': dimensions_path,
            },
            'field_counts_deriver': {
                # connect to a characteristic volume, which will remain constant
                'global': {
                    '_path': boundary_path,
                    'volume': ('characteristic_volume',)
                },
                'counts': ('species',),
                'concentrations': boundary_path + ('external',)
            },
            'delta_counts_to_concs': {
                'global': boundary_path,
                'counts': ('species',),
                'concentrations': ('delta_concentrations',),
            }
        }

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

        return topology



# plotting config

plot_variables_list_stochastic = [
    ('species', GLUCOSE_EXTERNAL),
    ('species', LACTOSE_EXTERNAL),
    ('species', 'rna_M'),
    ('species', 'protein_betaGal'),
    ('species', 'protein_Lactose_Permease'),
    ('flux_bounds', 'EX_glc__D_e'),
    ('flux_bounds', 'EX_lac__D_e'),
    ('boundary', ('mass', 'femtogram')),
    ('boundary', ('volume', 'femtoliter')),
]


# tests

def test_bioscrape_cobra_stochastic(
        total_time=2000,
        external_volume=1e-12 * units.L,
):
    bioscrape_composer = BioscrapeCOBRAstochastic({
        'local_fields': {'bin_volume': external_volume},
    })

    # get initial state
    initial_state = bioscrape_composer.initial_state()
    initial_state['boundary']['external'] = {
        GLUCOSE_EXTERNAL: 1e0,
        LACTOSE_EXTERNAL: 1e0}
    initial_state['boundary']['characteristic_volume'] = 10 * units.fL

    # make the experiment
    bioscrape_composite = bioscrape_composer.generate()
    bioscrape_experiment = Experiment(
        dict(
            processes=bioscrape_composite['processes'],
            topology=bioscrape_composite['topology'],
            initial_state=initial_state,))

    bioscrape_experiment.update(total_time)
    timeseries = bioscrape_experiment.emitter.get_timeseries()
    return timeseries

def run_bioscrape_cobra_stochastic(
    total_time=2000,
    out_dir='out',
):
    output = test_bioscrape_cobra_stochastic(total_time=total_time)

    # plot output
    variables_plot_config = {
        'out_dir': out_dir, 'filename': 'variables',
        'row_height': 2, 'row_padding': 0.2, 'column_width': 10,
        'variables': plot_variables_list_stochastic}

    plot_variables(output, **variables_plot_config)
    plot_simulation_output(output,
                           out_dir=out_dir,
                           filename='simulation_output')


def test_bioscrape_cobra_stochastic_divide(
        total_time=3000,
        external_volume=1e-12 * units.L,
):
    agent_id = '1'
    outer_path = ('agents', agent_id,)
    divide_config = {
        'divide_on': True,
        'agent_id': agent_id,
        'agents_path': ('..', '..', 'agents',),
        'fields_path': ('..', '..', 'fields',),
        'dimensions_path': ('..', '..', 'dimensions',),
        'local_fields': {'bin_volume': external_volume},
    }

    bioscrape_composer = BioscrapeCOBRAstochastic(divide_config)

    # get initial state
    initial_state = bioscrape_composer.initial_state()
    initial_state['boundary']['external'] = {
        GLUCOSE_EXTERNAL: 1e-1,
        LACTOSE_EXTERNAL: 1e-1}
    initial_state['boundary']['characteristic_volume'] = 10 * units.fL
    initial_state = {
        'agents': {
            agent_id: initial_state}}

    # make the experiment
    bioscrape_composite = bioscrape_composer.generate(path=outer_path)
    bioscrape_experiment = Experiment(
        dict(
            processes=bioscrape_composite['processes'],
            topology=bioscrape_composite['topology'],
            initial_state=initial_state,))

    bioscrape_experiment.update(total_time)
    timeseries = bioscrape_experiment.emitter.get_data_unitless()
    return timeseries

def run_bioscrape_cobra_stochastic_division(
        total_time=3000,
        out_dir='out'
):
    output = test_bioscrape_cobra_stochastic_divide(
        total_time=total_time)

    # multigen plots
    plot_settings = {
        'skip_paths': [
            # ('external',),
            ('internal_counts',),
            ('cobra_external',),
        ],
        'remove_zeros': False}
    plot_agents_multigen(
        output, plot_settings, out_dir, 'division_multigen')


# spatial test config
agent_id = '1'
outer_path = ('agents', agent_id,)
spatial_config = {
    'divide_on': True,
    'fields_on': True,
    'agent_id': agent_id,
    'agents_path': ('..', '..', 'agents',),
    'fields_path': ('..', '..', 'fields',),
    'dimensions_path': ('..', '..', 'dimensions',)}

# lattice environment test config
BOUNDS = [10, 10]
NBINS = [5, 5]
DEPTH = 10

def test_bioscrape_cobra_lattice(total_time=2500):

    # initial external
    field_concentrations = {
        GLUCOSE_EXTERNAL: 10,
        LACTOSE_EXTERNAL: 10,
    }

    # get initial state
    fields_composer = BioscrapeCOBRAstochastic(spatial_config)
    initial_state = fields_composer.initial_state()
    initial_state['boundary']['external'] = {
        GLUCOSE_EXTERNAL: 1e-1,
        LACTOSE_EXTERNAL: 1e-1}

    # initial agents
    initial_state = {
        'agents': {
            agent_id: initial_state}}

    # configure lattice compartment
    lattice_config_kwargs = {
        'bounds': BOUNDS,
        'n_bins': NBINS,
        'depth': DEPTH,
        'concentrations': field_concentrations}

    lattice_config = make_lattice_config(**lattice_config_kwargs)

    # declare the hierarchy
    hierarchy = {
        COMPOSER_KEY: {
            'type': Lattice,
            'config': lattice_config},
        'agents': {
            agent_id: {
                COMPOSER_KEY: {
                    'type': BioscrapeCOBRAstochastic,
                    'config': spatial_config}
            }}}

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


def run_bioscrape_cobra_stochastic_lattice(
        total_time=3000,
        out_dir='out'
):
    output = test_bioscrape_cobra_lattice(
        total_time=total_time
    )

    # multigen plots
    plot_settings = {
        'skip_paths': [
            # ('external',),
            ('internal_counts',),
            ('cobra_external',),
        ],
        'remove_zeros': True}
    plot_agents_multigen(
        output, plot_settings, out_dir, 'spatial_multigen')

    agents, fields = format_snapshot_data(output)
    plot_snapshots(
        bounds=BOUNDS,
        agents=agents,
        fields=fields,
        # include_fields=['glc__D_e', 'lcts_e'],
        include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
        out_dir=out_dir,
        filename='spatial_snapshots')

    tags_data = {
        'agents': agents,
        'fields': fields,
        'config': {'bounds': BOUNDS}}
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


def main():
    out_dir = os.path.join(EXPERIMENT_OUT_DIR, 'bioscrape_cobra_stochastic')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='bioscrape_cobra')
    parser.add_argument('--single', '-s', action='store_true', default=False)
    parser.add_argument('--divide', '-d', action='store_true', default=False)
    parser.add_argument('--fields', '-f', action='store_true', default=False)
    args = parser.parse_args()

    if args.single:
        biocobra_out_dir = os.path.join(out_dir, 'single')
        run_bioscrape_cobra_stochastic(
            total_time=2000,
            out_dir=biocobra_out_dir)

    if args.divide:
        div_out_dir = os.path.join(out_dir, 'division')
        run_bioscrape_cobra_stochastic_division(
            total_time=600,
            out_dir=div_out_dir)

    if args.fields:
        field_out_dir = os.path.join(out_dir, 'field')
        run_bioscrape_cobra_stochastic_lattice(
            total_time=600,
            out_dir=field_out_dir)

if __name__ == '__main__':
    main()
