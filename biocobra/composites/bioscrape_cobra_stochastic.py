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
from vivarium.core.composition import EXPERIMENT_OUT_DIR

# vivarium-bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# vivarium-cobra imports
from vivarium_cobra import Volume, LocalField
from vivarium_cobra.processes.configurations import get_iAF1260b_config
from vivarium_cobra.processes.dynamic_fba import DynamicFBA

# local imports
from biocobra.processes.flux_adaptor import DilutionFluxAdaptor, FluxAdaptor, AverageFluxAdaptor

# plots
from vivarium.plots.simulation_output import plot_simulation_output, plot_variables
from vivarium.plots.agents_multigen import plot_agents_multigen

NAME = 'BioscrapeCOBRA'
GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
SBML_FILE_DETERMINISTIC = 'lac_operon/LacOperon_deterministic.xml'
SBML_FILE_STOCHASTIC = 'lac_operon/LacOperon_stochastic.xml'

# choose the SBML file and set other bioscrape parameters
# deterministic_bioscrape_config = {
#     'sbml_file': SBML_FILE_DETERMINISTIC,
#     'stochastic': False,
#     'initial_volume': 1,
#     'internal_dt': 0.01,
# }
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
        'Lactose_consumed': {'input_type': 'delta'},  # No options specified
        'Glucose_internal': {'input_type': 'delta'},  # No options specified
    },
}
dilution_rate_flux_config = {
    'flux_keys': {
        'biomass': {
            'input_type': 'amount'
        }
    }
}
mass_mw_config = {
    'molecular_weights': {
        'mass': 1.0 * units.fg / units.molec
    }
}
# # convert stochastic delta counts to delta concentrations for constraining FBA
# delta_counts_to_concs_config = {
#     'keys': [
#         'Glucose_internal',
#         'Lactose_consumed',
#     ]
# }
# # configures the counts deriver to convert field concentrations to counts for stochastic sims
# field_counts_deriver_config = {
#     'keys': [GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL]}


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

# Here we override the default ports schema of the Biomass species and the k_dilution rate in Bioscrape.
# This is done so they can be set by the Derivers connected to mass and mass flux from Cobra.
schema_override = {
    'bioscrape': {
        'species': {
            'Biomass': {
                '_updater': 'set'  # override bioscrape ('species', 'Biomass') with a 'set' updater
            }
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
        # 'bioscrape_deterministic': deterministic_bioscrape_config,
        'bioscrape_stochastic': stochastic_bioscrape_config,
        'cobra': cobra_config,
        'flux_adaptor': flux_config,
        'dilution_rate_flux': dilution_rate_flux_config,
        # 'mass_to_molar': mass_mw_config,
        'mass_to_counts': mass_mw_config,
        'strip_units': strip_units_config,
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
        'stochastic': True,  # Is the CRN stochastic or deterministic?
        'spatial_on': False,  # are spatial dynamics used?
        'bioscrape_timestep': 1,
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
            'strip_units': StripUnits(config['strip_units'])}

        # Process Logic for different kinds of simulations

        # Stochastic Case
        # create a stochastic bioscrape model
        processes['bioscrape'] = Bioscrape(config['bioscrape_stochastic'])

        # flux is computed as an average flux
        processes['flux_adaptor'] = AverageFluxAdaptor(config['flux_adaptor'])

        # biomass is converted to a molecular count
        processes['biomass_adaptor'] = MassToCount(config['mass_to_counts'])

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
                'rates': {
                    '_path': ('rates',),
                },
                'globals': unitless_boundary_path,
            },
            'cobra': {
                'internal_counts': ('internal_counts',),
                # 'external': boundary_path + ('external',),
                'external': ('cobra_external',),  # These are handled separately from the external fields
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

        # # Ports added only in the stochastic case
        # else:
        #     pass

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

def test_bioscrape_cobra_stochastic(
        total_time=1000,
        external_volume=1e-12 * units.L,
):
    bioscrape_composer = BioscrapeCOBRAstochastic({
        'local_fields': {'bin_volume': external_volume},
    })

    # get initial state
    initial_state = bioscrape_composer.initial_state()
    initial_state['boundary']['external'] = {
        GLUCOSE_EXTERNAL: 1e2,
        LACTOSE_EXTERNAL: 1e2}

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

def test_bioscrape_cobra_stochastic_divide(
        total_time=1000,
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
        GLUCOSE_EXTERNAL: 1e2,
        LACTOSE_EXTERNAL: 1e2}
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
    timeseries = bioscrape_experiment.emitter.get_timeseries()
    return timeseries

# execute from the terminal
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

def run_bioscrape_cobra_stochastic(
    total_time=1000,
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
        'remove_zeros': True}
    plot_agents_multigen(
        output, plot_settings, out_dir, 'division_multigen')


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
            out_dir=div_out_dir)

    if args.fields:
        pass

if __name__ == '__main__':
    main()
