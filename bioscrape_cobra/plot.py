import matplotlib.pyplot as plt

# plotting
from vivarium.plots.simulation_output import plot_variables
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import (
    format_snapshot_data, plot_snapshots)
from vivarium_multibody.plots.snapshots import plot_tags


def plot_single(
        output,
        variables=None,
        out_dir=None,
        filename=None
):
    variables = variables or []
    variables_plot_config = {
        'out_dir': out_dir, 'filename': filename,
        'row_height': 2, 'row_padding': 0.2, 'column_width': 10,
        'variables': variables}
    fig = plot_variables(
        output, **variables_plot_config)
    return fig

def plot_multigen(
        output,
        variables=None,
        out_dir=None,
        filename=None,
):
    if filename and not out_dir:
        out_dir = 'out'

    plot_settings = {
        'include_paths': variables,
        'skip_paths': [
            ('internal_counts',),
            ('cobra_external',),
        ],
        'remove_zeros': False,
        'column_width': 8,
        'row_height': 1.5,
        'title_on_y_axis': True,
        'stack_column': True,
        'tick_label_size': 10,
        'title_size': 10}

    fig = plot_agents_multigen(
        output,
        plot_settings,
        out_dir,
        filename)
    return fig


def plot_fields(
        output,
        bounds=None,
        include_fields=None,
        tagged_molecules=None,
        n_snapshots=5,
        out_dir=None,
        filename=None
):
    agents, fields = format_snapshot_data(output)
    fig1 = plot_snapshots(
        bounds=bounds,
        agents=agents,
        fields=fields,
        include_fields=include_fields,
        out_dir=out_dir,
        filename=filename)

    tags_data = {
        'agents': agents,
        'fields': fields,
        'config': {'bounds': bounds}}
    tags_config = {
        'tagged_molecules': tagged_molecules,
        'n_snapshots': n_snapshots,
        'out_dir': out_dir,
        'filename': ('tags_' + filename) if filename else None}
    fig2 = plot_tags(
        data=tags_data,
        plot_config=tags_config)
    return fig1, fig2


# plotting
tags_dict = {
    'glc__D_e': 'tab:cyan',
    'co2_e': 'tab:orange',
    'h2o_c': 'tab:cyan',
    'atp_c': 'tab:orange',
    'asp__L_c': 'tab:green'
}

def move_to_end(data, d):
    for key in d.keys():
        if key in data:
            data[key] = data.pop(key)
    return data


def plot_metabolism(data, tags=tags_dict):
    # initialize subplots
    n_rows = 3
    n_cols = 1
    fig = plt.figure(figsize=(n_cols * 8, n_rows * 2))
    grid = plt.GridSpec(n_rows, n_cols)

    time_vec = data['time']

    # external
    ax = fig.add_subplot(grid[0, 0])
    external = move_to_end(data['external'], tags)
    for mol_id, series in external.items():
        if sum(series) != 0.0:
            ax.plot(
                time_vec,
                series,
                label=mol_id if mol_id in tags else None,
                color=tags.get(mol_id, 'tab:gray'),
            )
    ax.set_title(f'external ({len(external.keys())} species)')
    ax.set_ylabel('concentrations (mmol)')
    ax.set_yscale('log')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    # internal
    ax = fig.add_subplot(grid[1, 0])
    internal = move_to_end(data['internal_counts'], tags)
    for mol_id, series in internal.items():
        if sum(series) != 0.0:
            ax.plot(
                time_vec,
                series,
                label=mol_id if mol_id in tags else None,
                color=tags.get(mol_id, 'tab:gray'),
            )
    ax.set_title(f'internal ({len(internal.keys())} species)')
    ax.set_ylabel('counts')
    ax.set_yscale('log')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    # mass
    ax = fig.add_subplot(grid[2, 0])
    ax.plot(
        time_vec,
        data['global'][('mass', 'femtogram')],
        color='tab:blue',
    )
    ax.set_title('global')
    ax.set_ylabel('total mass (fg)')
    ax.set_xlabel('time (sec)')
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    fig.tight_layout()
    return fig
