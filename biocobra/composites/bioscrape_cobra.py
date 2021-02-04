import os
import numpy as np

# vivarium imports
from vivarium import TreeMass, DivideCondition, MetaDivision
from vivarium.core.process import Composite
from vivarium.library.units import units
from vivarium.core.composition import COMPOSITE_OUT_DIR, FACTORY_KEY, compose_experiment
from vivarium.core.experiment import Experiment

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
from biocobra.processes.flux_deriver import FluxDeriver
from biocobra.processes.biomass_adaptor import BiomassAdaptor

# plots
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data,
    plot_snapshots,
)
from vivarium_multibody.plots.snapshots import plot_tags


NAME = 'BioscrapeCOBRA'


#choose the SBML file and set other bioscrape parameters
deterministic_bioscrape_config = {
            'sbml_file': 'biocobra/data/LacOperon_deterministic.xml',
            'stochastic': False,
            'initial_volume': 1,
            'internal_dt': 0.01,}

# set cobra constrained reactions config
cobra_config = get_iAF1260b_config()
cobra_config.update({'time_step': 10})

#set up the config for the FluxDeriver
flux_config = {
    'flux_keys': [ #there are three fluxes keys, one for each flux computed
        'Lactose_consumed', 'Glucose_internal', 'biomass',
        ],
    'flux_options' : { #many options are toggled to get the correct biomass --> diluton rate conversion
        'biomass': [
            'amount', 'percent', 'positive', 'ignore zeros',
        ]
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
            'k_dilution__':{
                '_emit':True, #k_dilution should be emitted so it can be plotted
                '_updater':'set'
            }
        }
    }
}


class BioscrapeCOBRA(Composite):
    name = NAME
    defaults = {
        'bioscrape': deterministic_bioscrape_config,
        'cobra': cobra_config,
        'flux_deriver': flux_config,
        'divide_on': False,  # is division turned on?
        'agent_id': np.random.randint(0, 100),
        'divide_condition': {
            'threshold': 2000 * units.fg},
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'fields_path': ('..', '..', 'fields',),
        'dimensions_path': ('..', '..', 'dimensions',),
        'daughter_path': tuple(),
        '_schema': schema_override,
        'stochastic': False,
    }

    def generate_processes(self, config):

        bioscrape_process = Bioscrape(config['bioscrape'])
        cobra_process = DynamicFBA(config['cobra'])
        flux_deriver = FluxDeriver(config['flux_deriver'])

        # logic for counts/concentrations if stochastic/deterministic
        if config['bioscrape']['stochastic']:
            pass

        processes = {
            'bioscrape': bioscrape_process,
            'cobra': cobra_process,
            'mass_deriver': TreeMass(),
            'volume_deriver': Volume(),
            'flux_deriver': flux_deriver,
            'biomass_adaptor': BiomassAdaptor(),
            'local_field': LocalField(),
        }
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
                    # TODO (Eran) connect external molecules with boundary_path + ('external',)
                },
                'delta_species': ('delta_species',),
                'rates': {
                    '_path': ('rates',),
                    # 'k_dilution__': ('..', 'flux_bounds', 'k_dilution__'),
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
            'flux_deriver': {
                'deltas': ('delta_species',),
                'amounts': boundary_path,
                # connect Bioscrape deltas 'Lactose_consumed' and 'Glucose_internal'
                # to COBRA flux bounds 'EX_lac__D_e' and 'EX_glc__D_e'

                'fluxes':
                    {
                        '_path': ('flux_bounds',),
                        'Lactose_consumed': ('EX_lac__D_e',),
                        'Glucose_internal': ('EX_glc__D_e',),
                        'biomass': ('k_dilution__',)
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
            'local_field': {
                'exchanges': boundary_path + ('exchange',),
                'location': boundary_path + ('location',),
                'fields': fields_path,
                'dimensions': dimensions_path,
            },
        }

        # Ports added only in the deterministic case
        if not config['stochastic']:
            # connect biomass flux to the dilution rate
            topology['flux_deriver']['fluxes'].update({'biomass': ('k_dilution__',)})
            topology['bioscrape']['rates'].update({'k_dilution__': ('..', 'flux_bounds', 'k_dilution__')})

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
        return topology

agent_id = '1'
outer_path = ('agents', agent_id,)
fields_config = {
    'divide_on': True,
    'agent_id': agent_id,
    'agents_path': ('..', '..', 'agents',),
    'fields_path': ('..', '..', 'fields',),
    'dimensions_path': ('..', '..', 'dimensions',)}

def test_bioscrape_cobra():

    fields_composer = BioscrapeCOBRA(fields_config)

    initial_state = fields_composer.initial_state()
    initial_state['external']['glc__D_e'] = 10
    initial_state['external']['lcts_e'] = 10
    initial_state = {
        'agents': {
            agent_id: initial_state}}

    # make experiment
    fields_composite = fields_composer.generate(path=outer_path)
    fields_experiment = Experiment(
        dict(
            processes=fields_composite['processes'],
            topology=fields_composite['topology'],
            initial_state=initial_state,))

    total_time = 1000
    fields_experiment.update(total_time)
    output = fields_experiment.emitter.get_data()
    final_agents = output[total_time]['agents'].keys()
    assert len(final_agents) > 1, 'bioscrapeCOBRA agent did not successfully divide'
    return output

BOUNDS = [20, 20]
NBINS = [10, 10]
DEPTH = 10
def test_bioscrape_cobra_lattice(total_time=100):

    # get initial state
    fields_composer = BioscrapeCOBRA(fields_config)
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
                    'config': {}
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

    # output = test_bioscrape_cobra()
    output = test_bioscrape_cobra_lattice()

    plot_settings = {
        'remove_zeros': True}
    plot_agents_multigen(
        output, plot_settings, out_dir, 'bioCOBRA_lattice_multigen')

    agents, fields = format_snapshot_data(output)
    plot_snapshots(
        bounds=BOUNDS,
        agents=agents,
        fields=fields,
        include_fields=['glc__D_e', 'lcts_e'],
        out_dir=out_dir,
        filename='bioCOBRA_lattice_snapshots')

    tags_data = {'agents': agents, 'fields': fields, 'config': {'bounds': BOUNDS}}
    tags_config = {
        'tagged_molecules': [
            ('species', 'protein_Lactose_Permease',),
        ],
        'out_dir': out_dir,
        'filename': 'bioCOBRA_tags'}
    plot_tags(
        data=tags_data,
        plot_config=tags_config
        )



if __name__ == '__main__':
    run_bioscrape_cobra()
