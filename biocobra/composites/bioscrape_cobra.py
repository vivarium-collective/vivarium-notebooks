import numpy as np

# vivarium-core processes
from vivarium import TreeMass, Clock, MassToConcentration, MassToCount, DivideCondition, MetaDivision
from vivarium.core.process import Composer
from vivarium.library.units import units

# vivarium-bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# vivarium-cobra imports
from vivarium_cobra import Volume, LocalField
from vivarium_cobra.processes.configurations import get_iAF1260b_config
from vivarium_cobra.processes.dynamic_fba import DynamicFBA

# local imports
from biocobra.processes.flux_adaptor import DilutionFluxAdaptor, FluxAdaptor, AverageFluxAdaptor

NAME = 'BioscrapeCOBRA'
GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
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

#set cobra constrained reactions config
cobra_config = get_iAF1260b_config()

#set up the config for the FluxAdaptor
flux_config = {
    'flux_keys': {
        'Lactose_consumed': {'input_type': 'delta'}, #No options specified
        'Glucose_internal': {'input_type': 'delta'},  #No options specified
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
                '_updater': 'set',
            }
        }
    }
}


class BioscrapeCOBRA(Composer):
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
            processes['biomass_adaptor'] = MassToConcentration()

        # Stochastic Case
        else:
            # create a stochastic bioscrape model
            processes['bioscrape'] = Bioscrape(config['bioscrape_stochastic'])

            # flux is computed as an average flux
            processes['flux_adaptor'] = AverageFluxAdaptor(config['flux_adaptor'])

            # biomass is converted to a molecular count
            processes['biomass_adaptor'] = MassToCount()

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

        topology = {
            'bioscrape': {
                # all species go to a species store on the base level,
                # except Biomass, which goes to the 'globals' store, with variable 'biomass'
                'species': {
                    '_path': ('species',),
                    'Biomass': ('..',) + boundary_path + ('biomass',),
                    GLUCOSE_EXTERNAL: ('..',) + boundary_path + ('external', GLUCOSE_EXTERNAL,),
                    LACTOSE_EXTERNAL: ('..',) + boundary_path + ('external', LACTOSE_EXTERNAL,),
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
            topology.update({
                'local_field': {
                    'exchanges': boundary_path + ('exchange',),
                    'location': boundary_path + ('location',),
                    'fields': fields_path,
                    'dimensions': dimensions_path,
            }})

        return topology
