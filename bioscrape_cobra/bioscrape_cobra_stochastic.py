"""
================================
Stochastic Bioscrape/COBRA model
================================
"""
import os
import numpy as np

# vivarium-core processes
from vivarium import (
    TreeMass, Clock, MassToCount,
    CountsToMolar, DivideCondition,
    MetaDivision, MolarToCounts, StripUnits)
from vivarium.core.composer import Composer
from vivarium.library.units import units

# vivarium-bioscrape imports
from vivarium_bioscrape.processes.bioscrape import Bioscrape

# vivarium-cobra imports
from vivarium_cobra import Volume, LocalField
from vivarium_cobra.processes.configurations import get_iAF1260b_config
from vivarium_cobra.processes.cobra_fba import COBRA_FBA

# local imports
from bioscrape_cobra.flux_adaptor import AverageFluxAdaptor


GLUCOSE_EXTERNAL = 'Glucose_external'
LACTOSE_EXTERNAL = 'Lactose_external'
SBML_FILE_STOCHASTIC = 'LacOperon_stochastic.xml'
COBRA_TIMESTEP = 60
BIOSCRAPE_TIMESTEP = 20

# set mass threshold for division
divide_config = {
    'threshold': 2000 * units.fg}

# choose the SBML file and set other bioscrape parameters
stochastic_bioscrape_config = {
    'stochastic': True,
    'safe_mode': False,
    'initial_volume': 1,
    'internal_dt': 0.1}

# set cobra constrained reactions config
cobra_config = get_iAF1260b_config()
cobra_config["external"] = {"glc__D_e":{"_divider": "set"}}
cobra_config["external"] = {"lac__D_e":{"_divider": "set"}}

# set up the config for the FluxAdaptor
flux_config = {

    'flux_keys': {
        'Lactose_consumed': {
            'input_type': 'delta',
            'window_size': int(COBRA_TIMESTEP / BIOSCRAPE_TIMESTEP)},
        'Glucose_internal': {
            'input_type': 'delta',
            'window_size': int(COBRA_TIMESTEP / BIOSCRAPE_TIMESTEP)}}}

mass_mw_config = {
    'molecular_weights': {
        'mass': 1.0 * units.fg / units.molec}}

# convert stochastic delta counts to delta concentrations for constraining FBA
delta_counts_to_concs_config = {
    'keys': [
        'Glucose_internal',
        'Lactose_consumed']}

# configures the counts deriver to convert field concentrations to counts for stochastic sims
field_counts_deriver_config = {
    'keys': [GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL]}

# configuration for strip units deriver, which converts and removes specified units
strip_units_config = {
    'keys': [
        'mass', 'volume', 'density', 'biomass'],
    'convert': {
        # 'biomass': units.mmolar,
        'mass': units.ug}}

# Here we override the default ports schema of the Biomass species and the k_dilution rate in Bioscrape.
# This is done so they can be set by the Derivers connected to mass and mass flux from Cobra.
schema_override = {
    'bioscrape': {
        'species': {
            'Biomass': {
                # '_default': 0.00166,
                '_updater': 'set'},  # override bioscrape's Biomass with a 'set' updater
            'Glucose_external': {
                '_divider': 'set',
                '_updater': 'null'},
            'Lactose_external': {
                '_divider': 'set',
                '_updater': 'null'},
            'dna_Lac_Operon': {
                '_divider': 'set'}},
    },
    'field_counts_deriver': {
        'counts': {
            'Glucose_external': {
                '_divider': 'set'},
            'Lactose_external': {
                '_divider': 'set'}}
    },
}



# The stochastic Bioscrape/COBRA composer
class BioscrapeCOBRAstochastic(Composer):
    defaults = {
        'bioscrape_timestep': BIOSCRAPE_TIMESTEP,
        'cobra_timestep': COBRA_TIMESTEP,
        'divide_on': False,  # is division turned on?
        'fields_on': False,  # are spatial dynamics used?
        'reuse_processes': True,  # reuse the same processes for all agents?
        'sbml_file': SBML_FILE_STOCHASTIC,

        # process configs
        'bioscrape': stochastic_bioscrape_config,
        'cobra': cobra_config,
        'flux_adaptor': flux_config,
        'mass_to_counts': mass_mw_config,
        'delta_counts_to_concs': delta_counts_to_concs_config,
        'field_counts_deriver': field_counts_deriver_config,
        'strip_units': strip_units_config,
        'clock': {},
        'local_fields': {},

        # division config
        'agent_id': np.random.randint(0, 100),
        'divide_condition': divide_config,
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
        self.config['clock']['time_step'] = min(
            self.config['cobra_timestep'], self.config['bioscrape_timestep'])

        # configure parallelization
        self.config['cobra']['_parallel'] = self.config.get('_parallel', False)

        # sbml file
        self.config['bioscrape']['sbml_file'] = self.config['sbml_file']

        # configure local fields
        if not self.config['fields_on']:
            self.config['local_fields'].update({'nonspatial': True})

        # no processes initialized
        self.processes_initialized = False

    def initialize_processes(self, config):
        # Processes
        self.cobra = COBRA_FBA(config['cobra'])
        self.bioscrape = Bioscrape(config['bioscrape'])
        self.clock = Clock(config['clock'])

        # Derivers
        self.mass_deriver = TreeMass()
        self.volume_deriver = Volume()
        self.delta_counts_to_concs = CountsToMolar(config['delta_counts_to_concs'])
        self.biomass_adaptor = MassToCount(config['mass_to_counts'])
        self.strip_units = StripUnits(config['strip_units'])
        self.local_field = LocalField(config['local_fields'])
        self.field_counts_deriver = MolarToCounts(config['field_counts_deriver'])

        # set processes_initialized
        self.processes_initialized = self.config['reuse_processes']

    def generate_processes(self, config):
        if not self.processes_initialized:
            self.initialize_processes(config)

        processes = {
            # Processes
            'cobra': self.cobra,
            'bioscrape': self.bioscrape,
            'clock': self.clock,
            'flux_adaptor': AverageFluxAdaptor(config['flux_adaptor']),  # has internal state

            # Derivers
            'mass_deriver': self.mass_deriver,
            'volume_deriver': self.volume_deriver,
            'delta_counts_to_concs': self.delta_counts_to_concs,
            'biomass_adaptor': self.biomass_adaptor,
            'strip_units': self.strip_units,
            'local_field': self.local_field,
            'field_counts_deriver': self.field_counts_deriver}

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
        exchanges_path = boundary_path + ('exchanges',)

        topology = {
            'cobra': {
                'internal_counts': ('internal_counts',),
                'external': ('cobra_external',),
                'exchanges': {
                    # connect glc__D_e and lac__D_e to boundary exchanges that update fields
                    '_path': ('hidden_exchanges',),
                    'glc__D_e': ('..',) + exchanges_path + (GLUCOSE_EXTERNAL,),
                    'lac__D_e': ('..',) + exchanges_path + (LACTOSE_EXTERNAL,)},
                'reactions': ('reactions',),
                'flux_bounds': ('flux_bounds',),
                'global': boundary_path,
            },
            'bioscrape': {
                # all species go to a species store on the base level,
                # except Biomass, which goes to the 'boundary' store
                'species': {
                    '_path': ('species',),
                    'Biomass': ('..',) + unitless_boundary_path + ('biomass',)},
                'delta_species': ('delta_species',),
                'rates': ('rates',),
                'globals': unitless_boundary_path,
            },
            # convert the delta counts to delta concentrations
            'delta_counts_to_concs': {
                'global': boundary_path,
                'counts': ('delta_species',),
                'concentrations': ('delta_concentrations',),
            },
            # convert delta concentrations over several steps to an average flux bounds
            'flux_adaptor': {
                'inputs': ('delta_concentrations',),
                # connect Bioscrape deltas 'Lactose_consumed' and 'Glucose_internal'
                # to COBRA flux bounds 'EX_lac__D_e' and 'EX_glc__D_e'
                'fluxes': {
                    '_path': ('flux_bounds',),
                    'Lactose_consumed': ('EX_lac__D_e',),
                    'Glucose_internal': ('EX_glc__D_e',)}
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
            # apply exchanges from COBRA to the fields, and update external state
            'local_field': {
                'exchanges': exchanges_path,
                'location': boundary_path + ('location',),
                'fields': fields_path,
                'dimensions': dimensions_path,
            },
            # convert external concentrations to external counts, for Bioscrape to read from
            'field_counts_deriver': {
                # if fields_on, connect to value saved under (boundary, external),
                # otherwise connect directly to the fields_path, where a single value will be held
                'concentrations': boundary_path + ('external',) if config['fields_on'] else fields_path,
                'counts': ('species',),
                'global': {
                     # connect to a fixed bin volume
                     '_path': boundary_path,
                     'volume': ('bin_volume',)},
            },
        }

        if config['divide_on']:
            topology.update({
                'divide_condition': {
                    'variable': boundary_path + ('mass',),
                    'divide': boundary_path + ('divide',)},
                'division': {
                    'global': boundary_path,
                    'agents': agents_path}})

        return topology
