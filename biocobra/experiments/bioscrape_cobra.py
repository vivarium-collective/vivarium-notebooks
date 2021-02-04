import pylab as plt
import os

# vivarium imports
from vivarium.core.control import Control

from vivarium.core.composition import (
    simulate_process_in_experiment,
    simulate_compartment_in_experiment,
    compartment_in_experiment
)

# imported processes
from vivarium_cobra.processes.dynamic_fba import (
    DynamicFBA, get_iAF1260b_config)

from vivarium_bioscrape.processes.bioscrape import Bioscrape

# plots
from vivarium.plots.simulation_output import plot_simulation_output
from vivarium.plots.topology import plot_topology
from vivarium.plots.agents_multigen import plot_agents_multigen


# This works for deterministic/smooth dynamics, stochastic requires something else.
# calculate weighted average back in time?
# use delta since this was last used. get average delta based on COBRA's last update, use COBRA's timestep?
from biocobra.composites.bioscrape_cobra import BioscrapeCOBRA


# experiments
def run_bioscrape(
        total_time=2500,
        time_step=1,
):
    # initialize Bioscrape process
    bioscrape_config = {
        'sbml_file': 'Lac Operon Model/LacOperon_simple.xml',
        'time_step': time_step}
    bioscrape_process = Bioscrape(bioscrape_config)

    # initial state
    initial_state = bioscrape_process.initial_state()

    # run simulation
    settings = {
        'total_time': total_time,
        'initial_state': initial_state,
        'display_info': False,
        'progress_bar': False}
    bioscrape_timeseries = simulate_process_in_experiment(bioscrape_process, settings)

    return bioscrape_timeseries


def run_cobra(
        total_time=500,
        time_step=1,
):
    # get the configuration for the iAF1260b BiGG model
    config = get_iAF1260b_config()

    # load it into DynamicFBA
    metabolism = DynamicFBA(config)

    # get the model's initial state
    initial_state = metabolism.initial_state({})

    # run simulation
    sim_settings = {
        'initial_state': initial_state,
        'total_time': total_time}
    cobra_timeseries = simulate_process_in_experiment(metabolism, sim_settings)

    return cobra_timeseries



def run_bioscrape_cobra(
        total_time=2000,
        time_step=1,
):
    composite = BioscrapeCOBRA({})

    plot_topology(composite, {}, out_dir='out/experiments/bioscrape_cobra/', filename='network')

    initial_state = composite.initial_state()
    # initial_state['globals']['mass'] = initial_state['globals']['mass'].to('fg').magnitude

    # FBA external state
    # TODO connect these with Bioscrape, but don't let Bioscrape update them.
    initial_state['external']['glc__D_e'] = 10
    initial_state['external']['lcts_e'] = 10

    # run simulation
    sim_settings = {
        'initial_state': initial_state,
        'total_time': total_time}
    # experiment = compartment_in_experiment(composite, sim_settings)
    # experiment.update(total_time)
    output = simulate_compartment_in_experiment(composite, sim_settings)

    return output



def run_bioscrape_cobra_divide(
        total_time=2500,
        time_step=1,
):
    agent_id = '1'
    divide_config = {
        'divide_on': True,
        'agent_id': agent_id,
        'agents_path': ('..', '..', 'agents',),
    }
    division_composite = BioscrapeCOBRA(divide_config)

    plot_topology(division_composite, {}, out_dir='out/experiments/divide/', filename='network_divide')

    initial_state = division_composite.initial_state()
    initial_state['external']['glc__D_e'] = 10
    initial_state['external']['lcts_e'] = 10
    initial_state = {
            'agents': {
                agent_id: initial_state}}

    # run simulation
    # simulate
    settings = {
        'outer_path': ('agents', agent_id,),
        'initial_state': initial_state,
        'experiment_id': 'division'}
    division_experiment = compartment_in_experiment(
        division_composite,
        settings=settings,
        initial_state=initial_state)

    division_experiment.update(total_time)
    division_output = division_experiment.emitter.get_data_unitless()
    ts = division_experiment.emitter.get_timeseries()

    return division_output



def run_bioscrape_cobra_fields(
        total_time=100,
        time_step=1,
):
    agent_id = '1'
    fields_config = {
        'divide_on': True,
        'agent_id': agent_id,
        'agents_path': ('..', '..', 'agents',),
        'fields_path': ('..', '..', 'fields',),
        'dimensions_path': ('..', '..', 'dimensions',),
    }
    fields_composite = BioscrapeCOBRA(fields_config)

    plot_topology(fields_composite, {}, out_dir='out/experiments/fields/', filename='network_fields')

    initial_state = fields_composite.initial_state()
    initial_state['external']['glc__D_e'] = 10
    initial_state['external']['lcts_e'] = 10
    initial_state = {
            'agents': {
                agent_id: initial_state}}

    import ipdb;
    ipdb.set_trace()

    # # run simulation
    # # simulate
    # settings = {
    #     'outer_path': ('agents', agent_id,),
    #     'initial_state': initial_state,
    #     'experiment_id': 'division'}
    # fields_experiment = compartment_in_experiment(
    #     fields_composite,
    #     settings=settings,
    #     initial_state=initial_state)
    #
    # fields_experiment.update(total_time)
    # fields_output = fields_experiment.emitter.get_data_unitless()
    # ts = fields_experiment.emitter.get_timeseries()


    return fields_output



def plots_1(data, config, out_dir='out'):
    plot_simulation_output(
        data,
        settings={},
        out_dir=out_dir,
    )

def plots_2(data, config, out_dir='out'):
    plot_settings = {
        'skip_paths': [('external',)]
        #     'include_paths': [
        #         ('globals', 'mass'),
        #         ('species', 'Glucose_external'),
        #         ('internal_counts', 'zn2_c'),
        #         ('flux_bounds', 'EX_glc__D_e'),
        #         ('flux_bounds', 'EX_lac__D_e'),
        #         ],
    }
    plot_agents_multigen(
        data,
        settings=plot_settings,
        out_dir=out_dir,
    )



# plotting function for metabolism output
def plot_metabolism(data, config, out_dir='out'):

    ncol = config.get('ncol', 2)

    original_fontsize = plt.rcParams['font.size']
    plt.rcParams.update({'font.size': 9})

    # initialize subplots
    n_rows = 2
    n_cols = 2
    fig = plt.figure(figsize=(n_cols * 7, n_rows * 3))
    grid = plt.GridSpec(n_rows, n_cols)

    time_vec = data['time']

    # mass
    ax = fig.add_subplot(grid[0, 0])
    ax.plot(time_vec, data['global'][('mass', 'femtogram')], label='mass')
    ax.set_title('total compartment mass (fg)')
    ax.set_xlabel('time (sec)')
    #     ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=ncol)

    # external
    ax = fig.add_subplot(grid[0, 1])
    for mol_id, series in data['external'].items():
        if sum(series) != 0.0:
            ax.plot(time_vec, series, label=mol_id)
    ax.set_title('external concentrations (log)')
    ax.set_yscale('log')
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=ncol)

    # internal
    ax = fig.add_subplot(grid[1, 1])
    for mol_id, series in data['internal_counts'].items():
        if sum(series) != 0.0:
            ax.plot(time_vec, series, label=mol_id)
    ax.set_title('internal molecule counts (log)')
    ax.set_xlabel('time (sec)')
    ax.set_yscale('log')
    fig.tight_layout()
    plt.rcParams.update({'font.size': original_fontsize})

    os.makedirs(out_dir, exist_ok=True)
    filename = 'metabolism'
    # save figure
    fig_path = os.path.join(out_dir, filename)
    # plt.subplots_adjust(wspace=column_width / 3, hspace=column_width / 3)
    plt.savefig(fig_path, bbox_inches='tight')





# make the libraries
experiments_library = {
    '1': run_bioscrape,
    '2': run_cobra,
    '3': run_bioscrape_cobra,
    '4': run_bioscrape_cobra_divide,
    '5': run_bioscrape_cobra_fields,

}

plots_library = {
    '1': {
        'plot': plots_1,
        'config': {},
    },
    '2': {
        'plot': plot_metabolism,
        'config': {},
    },
    '3': {
        'plot': plots_2,
        'config': {},
    },
}

workflow_library = {
    'bioscrape': {
        'name': 'bioscrape',
        'experiment': '1',
        'plots': ['1'],
    },
    'cobra': {
        'name': 'cobra',
        'experiment': '2',
        'plots': ['2'],
    },
    'bioscrape_cobra': {
        'name': 'bioscrape_cobra',
        'experiment': '3',
        'plots': ['1'],
    },
    'divide': {
        'name': 'divide',
        'experiment': '4',
        'plots': ['3'],
    },
    'fields': {
        'name': 'divide',
        'experiment': '5',
        'plots': ['3'],
    },
}

if __name__ == '__main__':
    Control(
        experiments=experiments_library,
        plots=plots_library,
        workflows=workflow_library,
        out_dir='out/experiments',
        )
