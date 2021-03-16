"""
===================================
Deterministic Bioscrape/COBRA model
===================================
"""
import os
import numpy as np

# vivarium-core processes
from vivarium import (
    TreeMass, Clock, MassToMolar,
    DivideCondition, MetaDivision, StripUnits)
from vivarium.core.process import Composer
from vivarium.library.units import units

# vivarium-bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# vivarium-cobra imports
from vivarium_cobra import Volume, LocalField
from vivarium_cobra.processes.configurations import get_iAF1260b_config
from vivarium_cobra.processes.dynamic_fba import DynamicFBA

# local imports
from bioscrape_cobra.flux_adaptor import (
    DilutionFluxAdaptor, FluxAdaptor)

GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
dirname = os.path.dirname(__file__)
SBML_FILE_DETERMINISTIC = os.path.join(dirname, 'LacOperon_deterministic.xml')
COBRA_TIMESTEP = 50
BIOSCRAPE_TIMESTEP = 10

# choose the SBML file and set other bioscrape parameters
deterministic_bioscrape_config = {
    'sbml_file': SBML_FILE_DETERMINISTIC,
    'stochastic': False,
    'initial_volume': 1,
    'internal_dt': 0.01}

# set cobra constrained reactions config
cobra_config = get_iAF1260b_config()
cobra_config["external"] = {"glc__D_e":{"_divider": "set"}}
cobra_config["external"] = {"lac__D_e":{"_divider": "set"}}

# set up the config for the FluxAdaptor
flux_config = {
    'time_step': COBRA_TIMESTEP,
    'flux_keys': {
        'Lactose_consumed': {'input_type': 'delta'},  # No options specified
        'Glucose_internal': {'input_type': 'delta'}}}

dilution_rate_flux_config = {
    'time_step': COBRA_TIMESTEP,
    'flux_keys': {
        'mass': {
            'input_type': 'amount'}}}

mass_mw_config = {
    'molecular_weights': {
        'mass': 1.0 * units.fg / units.molec}}

# configuration for strip units deriver, which converts and removes specified units
strip_units_config = {
    'keys': [
        'mass', 'volume', 'density', 'biomass'],
    'convert': {
        'biomass': units.mmolar,
        'mass': units.ug}}

# Here we override the default ports schema of the Biomass species and the k_dilution rate in Bioscrape.
# This is done so they can be set by the Derivers connected to mass and mass flux from Cobra.
schema_override = {
    'bioscrape': {
        'species': {
            'Biomass': {
                '_default': 0.00166,
                '_updater': 'set'},  # override bioscrape ('species', 'Biomass') with a 'set' updater
            'Glucose_external': {
                '_divider': 'set',
                '_updater': 'null'},
            'Lactose_external': {
                '_divider': 'set',
                '_updater': 'null'},
            'Glucose_internal': {
                '_divider': 'set'},
            'Lactose_consumed': {
                '_divider': 'set'},
        },
        'rates': {
            'k_dilution__': {
                '_emit': True,  # k_dilution should be emitted so it can be plotted
                '_updater': 'set'}}}}


# The deterministic Bioscrape/COBRA composer
class BioscrapeCOBRAdeterministic(Composer):
    defaults = {
        'bioscrape_timestep': BIOSCRAPE_TIMESTEP,
        'cobra_timestep': COBRA_TIMESTEP,
        'divide_on': False,  # is division turned on?
        'fields_on': False,  # are spatial dynamics used?
        '_parallel': False,  # Are multiple cores used?
        'reuse_processes': True,  # reuse the same processes for all agents?

        # process configs
        'bioscrape': deterministic_bioscrape_config,
        'cobra': cobra_config,
        'flux_adaptor': flux_config,
        'dilution_rate_flux': dilution_rate_flux_config,
        'mass_to_molar': mass_mw_config,
        'strip_units': strip_units_config,
        'clock': {},
        'local_fields': {
            'bin_volume': 1e-12 * units.L},

        # division config
        'agent_id': np.random.randint(0, 1000),
        'divide_condition': {
            'threshold': 2000 * units.fg},
        'boundary_path': ('boundary',),
        'agents_path': ('agents',),
        'fields_path': ('fields',),
        'dimensions_path': ('dimensions',),
        'daughter_path': tuple(),

        # override
        '_schema': schema_override}

    def __init__(self, config=None):
        super().__init__(config)

        # configure timesteps
        self.config['bioscrape']['time_step'] = self.config['bioscrape_timestep']
        self.config['flux_adaptor']['time_step'] = self.config['bioscrape_timestep']
        self.config['cobra']['time_step'] = self.config['cobra_timestep']
        self.config['dilution_rate_flux']['time_step'] = self.config['cobra_timestep']
        self.config['clock']['time_step'] = min(
            self.config['cobra_timestep'], self.config['bioscrape_timestep'])

        #configure parallelization
        self.config['cobra']['_parallel'] = self.config.get('_parallel', False)

        # configure local fields
        if not self.config['fields_on']:
            self.config['local_fields'].update({'nonspatial': True})

        # no processes initialized
        self.processes_initialized = False

    def initialize_processes(self, config):
        # Processes
        self.cobra = DynamicFBA(config['cobra'])
        self.bioscrape = Bioscrape(config['bioscrape'])
        self.clock = Clock(config['clock'])

        # Derivers
        self.mass_deriver = TreeMass()
        self.volume_deriver = Volume()
        self.biomass_adaptor = MassToMolar(config['mass_to_molar'])
        self.strip_units = StripUnits(config['strip_units'])
        self.local_field = LocalField(config['local_fields'])

        if self.config['reuse_processes']:
            self.processes_initialized = True

    def generate_processes(self, config):
        if not self.processes_initialized:
            self.initialize_processes(config)

        processes = {
            # Processes
            'cobra': self.cobra,
            'bioscrape': self.bioscrape,
            'flux_adaptor': FluxAdaptor(config['flux_adaptor']),  # has internal state
            'dilution_rate_adaptor': DilutionFluxAdaptor(config['dilution_rate_flux']),  # has internal state
            'clock': self.clock,

            # Derivers
            'mass_deriver': self.mass_deriver,
            'volume_deriver': self.volume_deriver,
            'biomass_adaptor': self.biomass_adaptor,
            'strip_units': self.strip_units,
            'local_field': self.local_field}

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
            'cobra': {
                'internal_counts': ('internal_counts',),
                'external': ('cobra_external',),
                'exchanges': {
                    # connect glc__D_e and lac__D_e to boundary exchanges that update fields
                    '_path': ('hidden_exchanges',),
                    'glc__D_e': ('..',) + boundary_path + ('exchanges', GLUCOSE_EXTERNAL,),
                    'lac__D_e': ('..',) + boundary_path + ('exchanges', LACTOSE_EXTERNAL,)},
                'reactions': ('reactions',),
                'flux_bounds': ('flux_bounds',),
                'global': boundary_path,
            },
            'bioscrape': {
                # all species go to a species store on the base level,
                # except Biomass, which goes to the 'boundary' store
                'species': {
                    '_path': ('species',),
                    'Biomass': ('..',) + unitless_boundary_path + ('biomass',),
                    GLUCOSE_EXTERNAL: ('..',) + boundary_path + ('external', GLUCOSE_EXTERNAL,),
                    LACTOSE_EXTERNAL: ('..',) + boundary_path + ('external', LACTOSE_EXTERNAL,)},
                'delta_species': ('delta_species',),
                'rates': ('rates',),
                'globals': unitless_boundary_path,
            },
            'local_field': {
                'exchanges': boundary_path + ('exchanges',),
                'location': boundary_path + ('location',),
                # connect fields directly to external port if fields_on is False
                'fields': fields_path if config['fields_on'] else boundary_path + ('external',),
                'dimensions': dimensions_path,
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
                    'mass': ('mass',)},
                'output': {
                    '_path': boundary_path,
                    'mass': ('biomass',)},
                'global': boundary_path,
            },
            'clock': {
                'global_time': boundary_path + ('time',)
            },
            'strip_units': {
                'units': boundary_path,
                'no_units': unitless_boundary_path,
            },
            # Create port biomass flux to the dilution rate computed by the dilution_rate_adaptor process
            'dilution_rate_adaptor': {
                'inputs': boundary_path,
                'fluxes': {
                    '_path': ('rates',),
                    'mass': ('k_dilution__',)}
            }
        }

        if config['divide_on']:
            topology.update({
                'divide_condition': {
                    'variable': boundary_path + ('mass',),
                    'divide': boundary_path + ('divide',),},
                'division': {
                    'global': boundary_path,
                    'agents': agents_path}})

        return topology
