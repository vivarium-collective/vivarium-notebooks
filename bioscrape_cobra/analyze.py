"""
Analyze an experiment saved on the database emitter (mongoDB) by running:
    `python bioscrape_cobra/analyze.py experiment_id`

These functions were specifically built to generate Figure 9, which is a stochastic, spatial simulation

To access the saved experiments on the mongoDB database, install vivarium scripts with:
    `pip install vivarium-scripts`
and run the commands:
    `python -m scripts.access_db list`
    `python -m scripts.access_db info [experiment_id]`
"""

import os
import argparse
import copy

from vivarium.core.emitter import (
    data_from_database,
    DatabaseEmitter,
    remove_units, deserialize_value)
from vivarium.plots.simulation_output import save_fig_to_dir, set_axes
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data, make_tags_figure, get_agent_colors)
from bioscrape_cobra.plot import (
    plot_fields_tags, plot_fields_snapshots)
from bioscrape_cobra.bioscrape_cobra_stochastic import (
    GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL)


MULTIGEN_PLOT_CONFIG = {
    'include_paths': [
        ('species', 'rna_M'),
        ('species', 'protein_betaGal'),
        ('species', 'protein_Lactose_Permease'),
        ('species', 'Glucose_external'),
        ('species', 'Lactose_external'),
        # ('flux_bounds', 'EX_glc__D_e'),
        # ('flux_bounds', 'EX_lac__D_e'),
        ('boundary', 'volume'),
        ('boundary', 'mass'),
    ],
    'store_order': ('species', 'flux_bounds', 'boundary'),
    'titles_map': {
        ('species', 'rna_M'): 'lac operon RNA',
        ('species', 'protein_betaGal'): r'$\beta$-galactosidase',
        ('species', 'protein_Lactose_Permease'): 'lactose permease',
        ('species', 'Glucose_external'): 'external glucose',
        ('species', 'Lactose_external'): 'external lactose',
        ('flux_bounds', 'EX_glc__D_e'): 'glucose flux bound (mmol/L/s)',
        ('flux_bounds', 'EX_lac__D_e'): 'lactose flux bound (mmol/L/s)',
        ('boundary', 'volume'): 'volume (fL)',
    },
    'remove_zeros': False,
    'column_width': 5,
    'row_height': 1.5,
    'title_on_y_axis': False,
    'stack_column': True,
    'tick_label_size': 10,
    'title_size': 10,
    'sci_notation': 3,
}

YLABEL_SIZE = 48

def access(experiment_id):

    # mongo client
    config = {
        'host': '{}:{}'.format('localhost', 27017),
        'database': 'simulations'}
    emitter = DatabaseEmitter(config)
    db = emitter.db

    # access experiment from database
    data, experiment_config = data_from_database(experiment_id, db)

    # get the bounds
    bounds = experiment_config['processes']['diffusion']['bounds']

    # reformat data
    deserialized = deserialize_value(data)
    output = remove_units(deserialized)

    return output, bounds


def plot_fields_fig(output, bounds, out_dir):
    new_output = copy.deepcopy(output)
    seconds_list = list(new_output.keys())
    for seconds in seconds_list:
        hours = round(seconds/60/60, 2)
        new_output[hours] = new_output.pop(seconds)

    # plot snapshots fields
    fig_snapshots = plot_fields_snapshots(
        new_output,
        bounds=bounds,
        include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
        colorbar_decimals=1,
        agent_fill_color='k',
        show_timeline=True,
        time_unit='hr')

    # alter figure labels and save
    axes = fig_snapshots.get_axes()
    for axis in axes:
        axis.set_title('')
        ylabel = axis.get_ylabel()
        if ylabel == 'Glucose_external':
            axis.set_ylabel('')
            axis.set_title('external glucose field', fontsize=YLABEL_SIZE, pad=15, x=2.3)
        if ylabel == 'Lactose_external':
            axis.set_ylabel('')
            axis.set_title('external lactose field', fontsize=YLABEL_SIZE, pad=15, x=2.3)

    save_fig_to_dir(
        fig_snapshots,
        out_dir=out_dir,
        filename='field_snapshots.pdf')


def plot_phylogeny_fig(output, bounds, agent_colors=None, out_dir='out'):

    # plot phylogeny snapshots
    fig_phylogeny = plot_fields_snapshots(
        output,
        bounds=bounds,
        skip_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
        agent_colors=agent_colors,
        colorbar_decimals=1,
        show_timeline=False,
        filename='phylogeny_snapshots.pdf',
        out_dir=out_dir)


def plot_tags_fig(output, bounds, out_dir):
    # plot tags
    fig_tags = plot_fields_tags(
        output,
        bounds=bounds,
        convert_to_concs=False,
        tagged_molecules=[('species', 'protein_Lactose_Permease',)],
        colorbar_decimals=1)
    # alter figure and save
    axes = fig_tags.get_axes()
    for axis in axes:
        axis.set_title('')
        ylabel = axis.get_ylabel()
        if ylabel:
            axis.set_ylabel('internal\nlactose\npermease', fontsize=YLABEL_SIZE)
    save_fig_to_dir(
        fig_tags,
        out_dir=out_dir,
        filename='tags_snapshots.pdf')


def plot_multigen_fig(output, agent_colors=None, out_dir='out'):
    plot_config = copy.deepcopy(MULTIGEN_PLOT_CONFIG)
    plot_config['agent_colors'] = agent_colors

    # plot multigen
    multigen_fig = plot_agents_multigen(output, plot_config)

    # add colony mass to multigen_fig
    colony_mass, time_vec = get_colony_mass(output)

    # convert to hours
    time_vec = [time/3600 for time in time_vec]

    # get mass axis to replace with colony mass
    allaxes = multigen_fig.get_axes()
    ax = None
    for axis in allaxes:
        if axis.get_title() == 'boundary \nmass':
            ax = axis

    # replace with colony mass
    ax.clear()
    set_axes(ax, True, sci_notation=MULTIGEN_PLOT_CONFIG['sci_notation'])
    ax.plot(time_vec, colony_mass, linewidth=3.0, color='darkslategray')
    ax.set_xlim([time_vec[0], time_vec[-1]])
    ax.set_title('total colony mass (fg)', rotation=0, fontsize=MULTIGEN_PLOT_CONFIG['title_size'])
    ax.set_xlabel('time (hr)')
    ax.spines['bottom'].set_position(('axes', -0.2))

    # save
    save_fig_to_dir(multigen_fig, 'spatial_multigen.pdf', out_dir)


def plot_single_tags(agents, bounds, out_dir):

    # make individual tag plots
    time_vec = list(agents.keys())
    snapshot_times = [time_vec[-1]]
    time_indices = [time_vec.index(time) for time in snapshot_times]

    ############
    # glc flux #
    ############
    fig_tags = make_tags_figure(
        agents=agents, bounds=bounds, n_snapshots=1,
        time_indices=time_indices, snapshot_times=snapshot_times,
        background_color='white', scale_bar_length=False, show_timeline=False, convert_to_concs=False,
        tagged_molecules=[('flux_bounds', 'EX_glc__D_e')])

    # alter figure and save
    axes = fig_tags.get_axes()
    axes[0].set_title('glucose flux\n(mmol/L/s)', fontsize=YLABEL_SIZE, pad=15)
    axes[0].get_yaxis().set_visible(False)
    save_fig_to_dir(fig_tags, out_dir=out_dir, filename='tag_EX_glc__D_e.pdf')


    ############
    # lac flux #
    ############
    fig_tags = make_tags_figure(
        agents=agents, bounds=bounds, n_snapshots=1,
        time_indices=time_indices, snapshot_times=snapshot_times,
        background_color='white', scale_bar_length=False, show_timeline=False, convert_to_concs=False,
        tagged_molecules=[('flux_bounds', 'EX_lac__D_e')])
    # alter figure and save
    axes = fig_tags.get_axes()
    axes[0].set_title('lactose flux\n(mmol/L/s)', fontsize=YLABEL_SIZE, pad=15)
    axes[0].get_yaxis().set_visible(False)
    save_fig_to_dir(fig_tags, out_dir=out_dir, filename='tag_EX_lac__D_e.pdf')


    ################
    # lac Permease #
    ################
    fig_tags = make_tags_figure(
        agents=agents, bounds=bounds, n_snapshots=1,
        time_indices=time_indices, snapshot_times=snapshot_times,
        background_color='white', scale_bar_length=False, show_timeline=False, convert_to_concs=False,
        tagged_molecules=[('species', 'protein_Lactose_Permease')])
    # alter figure and save
    axes = fig_tags.get_axes()
    axes[0].set_title('lactose permease\n(counts)', fontsize=YLABEL_SIZE, pad=15)
    axes[0].get_yaxis().set_visible(False)
    save_fig_to_dir(fig_tags, out_dir=out_dir, filename='tag_protein_Lactose_Permease.pdf')


    ###############
    # growth rate #
    ###############
    agents = add_growth_rate_to_agents(agents)
    fig_tags = make_tags_figure(
        agents=agents, bounds=bounds, n_snapshots=1,
        time_indices=time_indices, snapshot_times=snapshot_times,
        background_color='white', scale_bar_length=False, show_timeline=False, convert_to_concs=False,
        tagged_molecules=[('boundary', 'growth_rate')])
    # alter figure and save
    axes = fig_tags.get_axes()
    axes[0].set_title('growth rate\n(fg/s)', fontsize=YLABEL_SIZE, pad=15)
    axes[0].get_yaxis().set_visible(False)
    save_fig_to_dir(fig_tags, out_dir=out_dir, filename='tag_growth_rate.pdf')


def add_growth_rate_to_agents(agents):
    time_vec = list(agents.keys())

    # initial state at 0
    for agent_id, state in agents[time_vec[0]].items():
        state['boundary']['growth_rate'] = 0

    for t0, t1 in zip(time_vec[:-1], time_vec[1:]):
        agents0 = agents[t0]
        agents1 = agents[t1]
        for agent_id, state1 in agents1.items():
            if agent_id in agents0:
                mass0 = agents0[agent_id]['boundary']['mass']
                mass1 = state1['boundary']['mass']
                dm = mass1 - mass0
                dt = t1 - t0
                growth_rate = dm/dt
                state1['boundary']['growth_rate'] = growth_rate
    return agents


def get_colony_mass(output):
    colony_mass = []
    time_vec = list(output.keys())
    for time, state in output.items():
        mass = 0.0
        for agent_id, agent_state in state['agents'].items():
            mass += agent_state['boundary']['mass']
        colony_mass.append(mass)
    return colony_mass, time_vec


def main():

    # parse
    parser = argparse.ArgumentParser(description='access data from db')
    parser.add_argument('experiment_id', type=str, default=False)
    parser.add_argument('--multigen', '-1', action='store_true', default=False)
    parser.add_argument('--fields', '-2', action='store_true', default=False)
    parser.add_argument('--tags', '-3', action='store_true', default=False)
    parser.add_argument('--single_tags', '-4', action='store_true', default=False)
    parser.add_argument('--all', '-a', action='store_true', default=False)
    args = parser.parse_args()
    experiment_id = args.experiment_id


    # make a directory for the figures
    out_dir = f'out/analyze/{experiment_id}'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # retrieve the data
    output, bounds = access(experiment_id)
    del output[0.0]

    agents, fields = format_snapshot_data(output)

    # run the plot functions
    if args.multigen or args.all:
        agent_colors = get_agent_colors(agents)
        plot_phylogeny_fig(output, bounds, agent_colors, out_dir)
        plot_multigen_fig(output, agent_colors, out_dir)

    if args.fields or args.all:
        plot_fields_fig(output, bounds, out_dir)

    if args.tags or args.all:
        plot_tags_fig(output, bounds, out_dir)

    if args.single_tags or args.all:
        plot_single_tags(agents, bounds, out_dir)

if __name__ == '__main__':
    main()
