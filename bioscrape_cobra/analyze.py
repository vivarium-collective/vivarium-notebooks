"""
Analyze an experiment saved on the database emitter (mongoDB) by running:
    `python bioscrape_cobra/analyze.py experiment_id`

To access the saved experiments on the mongoDB database, install vivarium scripts with:
    `pip install vivarium-scripts`
and run the commands:
    `python -m scripts.access_db list`
    `python -m scripts.access_db info [experiment_id]`
"""

import os
import argparse

import matplotlib.pyplot as plt

from vivarium.core.emitter import (
    data_from_database,
    DatabaseEmitter,
    remove_units, deserialize_value
)
from vivarium.plots.simulation_output import save_fig_to_dir
from bioscrape_cobra.plot import (plot_multigen, plot_single, plot_fields_tags, plot_fields_snapshots)
from bioscrape_cobra.bioscrape_cobra_stochastic import (
    GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL)

MULTIGEN_VARIABLES = [
    ('species', 'rna_M'),
    ('species', 'protein_betaGal'),
    ('species', 'protein_Lactose_Permease'),
    ('flux_bounds', 'EX_glc__D_e'),
    ('flux_bounds', 'EX_lac__D_e'),
    ('boundary', 'volume'),
]

MULTIGEN_PLOT_CONFIG = {
    'store_order': ('species', 'flux_bounds', 'boundary'),
    'titles_map': {
        ('species', 'rna_M'): 'rna\nM',
        ('species', 'protein_betaGal'): 'betaGal',
        ('species', 'protein_Lactose_Permease'): 'Lactose\nPermease',
        ('flux_bounds', 'EX_glc__D_e'): 'glc\nflux\nbound',
        ('flux_bounds', 'EX_lac__D_e'): 'lac\nflux\nbound',
        ('boundary', 'volume'): 'volume',
    },
}

def access():

    # parse
    parser = argparse.ArgumentParser(description='access data from db')
    parser.add_argument(
            'experiment_id',
            type=str,
            default=False)
    args = parser.parse_args()
    experiment_id = args.experiment_id

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

    return output, bounds, experiment_id


def plot_full(output, bounds, experiment_id):
    # make a directory for the figures
    out_dir = f'out/analyze/{experiment_id}'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # plot
    multigen_fig = plot_multigen(
        output,
        variables=MULTIGEN_VARIABLES,
        out_dir=out_dir,
        filename='spatial_multigen',
        **MULTIGEN_PLOT_CONFIG)

    fig_snapshots = plot_fields_snapshots(
        output,
        bounds=bounds,
        include_fields=[GLUCOSE_EXTERNAL, LACTOSE_EXTERNAL],
        colorbar_decimals=1,
        agent_fill_color='k',
        show_timeline=False,
    )

    ylabel_size = 48
    # alter figure and save
    axes = fig_snapshots.get_axes()
    for axis in axes:
        axis.set_title('')
        ylabel = axis.get_ylabel()
        if ylabel == 'Glucose_external':
            axis.set_ylabel('external\nglucose', fontsize=ylabel_size)
        if ylabel == 'Lactose_external':
            axis.set_ylabel('external\nlactose', fontsize=ylabel_size)

    size = fig_snapshots.get_size_inches()
    save_fig_to_dir(
        fig_snapshots,
        out_dir=out_dir,
        filename='bioscrape_cobra_stochastic_lattice_snapshots.pdf')

    # tags figure
    fig_tags = plot_fields_tags(
        output,
        bounds=bounds,
        tagged_molecules=[('species', 'protein_Lactose_Permease',)],
        colorbar_decimals=1,
    )

    # alter figure and save
    axes = fig_tags.get_axes()
    for axis in axes:
        axis.set_title('')
        ylabel = axis.get_ylabel()
        if ylabel:
            axis.set_ylabel('internal\nlactose\npermease', fontsize=ylabel_size)

    fig_snapshots.set_size_inches(size[0], size[1]/2, forward=True)
    save_fig_to_dir(
        fig_tags,
        out_dir=out_dir,
        filename='bioscrape_cobra_stochastic_lattice_tags.pdf')



if __name__ == '__main__':
    output, bounds, experiment_id = access()
    plot_full(output, bounds, experiment_id)
