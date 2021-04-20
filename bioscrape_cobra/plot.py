import matplotlib.pyplot as plt

# plotting
from vivarium.plots.simulation_output import plot_variables
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_multibody.plots.snapshots import format_snapshot_data, plot_snapshots, plot_tags
import copy

def plot_single(
        output,
        variables=None,
        config=None,
        out_dir=None,
        filename=None
):
    variables = variables or []
    config = config or {}
    variables_plot_config = {
        'out_dir': out_dir,
        'filename': filename,
        'row_height': 1.2,
        'row_padding': 0.6,
        'column_width': 4,
        'default_color': 'tab:gray',
        'linewidth': 2.0,
        'sci_notation': 3,
        'variables': variables}
    variables_plot_config.update(config)
    fig = plot_variables(
        output, **variables_plot_config)
    return fig

def plot_multigen(
        output,
        variables=None,
        out_dir=None,
        filename=None,
        **kwargs
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
        'title_size': 10,
        **kwargs
        }

    fig = plot_agents_multigen(
        output,
        plot_settings,
        out_dir,
        filename)
    return fig


def plot_fields_snapshots(
        output,
        bounds=None,
        include_fields=None,
        n_snapshots=4,
        colorbar_decimals=4,
        scale_bar_length=5,
        phylogeny_colors=False,
        out_dir=None,
        filename=None,
        **kwargs):

    agents, fields = format_snapshot_data(output)

    fig1 = plot_snapshots(
        bounds=bounds,
        agents=agents,
        fields=fields,
        phylogeny_names=phylogeny_colors,
        n_snapshots=n_snapshots,
        scale_bar_length=scale_bar_length,
        colorbar_decimals=colorbar_decimals,
        include_fields=include_fields,
        out_dir=out_dir,
        filename=filename,
        **kwargs)

    return fig1

# plotting
tags_dict = {
    'glc__D_e': 'tab:cyan',
    'co2_e': 'tab:orange',
    'h2o_c': 'tab:cyan',
    'atp_c': 'tab:orange',
    'asp__L_c': 'tab:green'
}
def plot_fields_tags(
        output,
        bounds=None,
        tagged_molecules=None,
        n_snapshots=4,
        colorbar_decimals=4,
        out_dir=None,
        filename=None,
        **kwargs):

    fig2 = plot_tags(
        data=output,
        bounds=bounds,
        tagged_molecules=tagged_molecules,
        n_snapshots=n_snapshots,
        background_color='w',
        scale_bar_length=5,
        colorbar_decimals=colorbar_decimals,
        out_dir=out_dir,
        filename=('tags_' + filename) if filename else None,
        **kwargs
    )
    
    return fig2




def move_to_end(data, d):
    for key in d.keys():
        if key in data:
            data[key] = data.pop(key)
    return data


def plot_metabolism(data, tags=tags_dict):
    # initialize subplots

    n_rows = 3
    n_cols = 1
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 4/3))
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
        color='tab:blue', label="biomass (fg)"
    )
    ax.set_ylabel('cell mass')
    ax.set_title('global')
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    #Set the x label for the last axis, whichever it is
    ax.set_xlabel('time (sec)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    fig.tight_layout()
    return fig


# Below are different configs for the topology plots
custom_widths = {
    'cobra':8,
    'bioscrape':8,
    'division':8,
    'diffusion':8,
    'multibody':8
}

embedded_custom_widths = {"agents\n0\n"+k:custom_widths[k] for k in custom_widths if k not in ["diffusion", "multibody"]}
embedded_custom_widths['diffusion'] = custom_widths['diffusion']
embedded_custom_widths['multibody'] = custom_widths['multibody']

u = 1.0
r0 = 0
r1 = 1.5*u
cobra_topology_config = {
    'graph_format': 'hierarchy',
    'dashed_edges': True,
    'show_ports': False,
    'remove_nodes':[
        'exchanges', 'field_deriver', 'global\nlocation', 'dimensions',
    ],
    'process_colors':{
        'cobra':'tab:orange',
        'mass\nderiver':'tab:orange',
        'volume\nderiver':'tab:orange',
    },
    'store_colors':{
        'internal_counts':'tab:orange',
        'external':'tab:orange',
        'exchanges':'tab:orange',
        'reactions':'tab:orange',
        'flux_bounds':'tab:orange',
        'global':'tab:brown',
    },
    'node_labels':{
        'mass_deriver':'mass\nderiver',
        'volume_deriver':'volume\nderiver',
        'internal_counts':'internal\ncounts',
        'flux_bounds':'flux\nbounds',
    },
    'coordinates':{
        
        'cobra':(-1*u, r1),
        'mass_deriver':(0*u, r1),
        'volume_deriver':(1*u, r1),
        
        'reactions':(-2*u, r0),
        'internal_counts':(-1*u, r0),
        'external':(0*u, r0),
        'flux_bounds':(1*u, r0),
        'global':(2*u, r0)
    },
    'custom_widths':custom_widths
}

bioscrape_topology_config = {
    'dashed_edges':True,
    'show_ports':False,
    'process_colors':{
        'bioscrape':'tab:green',
    },
    'store_colors':{
        'species':'tab:green',
        'delta_species':'tab:green',
        'rates':'tab:green',
        'globals':'tab:brown'
    },
    'node_labels':{
        'delta_species':'delta\nspecies'
    },
    'custom_widths':custom_widths
}


u = 0.9
config_grow_divide_topology_i = {
    'graph_format': 'hierarchy',
    'dashed_edges':True,
    'show_ports':False,
    'remove_nodes':[
        'agents\n0\nrates\ngrowth_rate', 'agents\n0\nrates\ngrowth_noise', 'dimensions'
    ],
    'process_colors':{
        'agents\n0\ndivide_condition':'tab:cyan',
        'agents\n0\ndivision':'tab:cyan',
        'agents\n0\nglobals_deriver':'tab:gray',
        'agents\n0\ngrowth':'tab:gray',
    },
    'store_colors':{
        'agents':'tab:cyan',
        'agents\n0':'tab:cyan',
        'agents\n0\nrates':'tab:gray',
        'agents\n0\nboundary':'tab:brown',
        'agents\n0\nboundary\nmass':'tab:brown',
        'agents\n0\nboundary\ndivide':'tab:brown',
    },
    'coordinates':{
        'agents\n0\ndivision':(u, -3*u),
        'agents\n0\ndivide_condition':(2*u, -3*u),
        
        'agents':(0,  u),
        'agents\n0':(0, 0),
        'agents\n0\nrates':(-1*u, -u),
        'agents\n0\nboundary':(0*u, -u),
        'agents\n0\nboundary\nmass':(1*u, -1.5*u),
        'agents\n0\nboundary\ndivide':(2*u, -1.5*u),
        
        'agents\n0\ngrowth':(-1*u, -3*u),
        'agents\n0\nglobals_deriver':(0*u, -3*u)
    },
    'node_labels':{
        'agents\n0':"agent\n0",
        'agents\n0\nglobals_deriver':'agent 0\nglobals\nderiver',
        'agents\n0\ndivide_condition':'agent 0\ndivide\ncondition',
        'agents\n0\ndivision':'agent 0\ndivision',
        'agents\n0\nrates':'agent 0\nrates',
        'agents\n0\ngrowth':'agent 0\ngrowth',
        'agents\n0\ngloabls\nderiver':'agent 0\nglobals\nderiver',
        'agents\n0\nboundary':'agent 0\nboundary',
        'agents\n0\nboundary\nmass':'agent 0\nboundary\nmass',
        'agents\n0\nboundary\ndivide':'agent 0\nboundary\ndivide',
    },
    'custom_widths':embedded_custom_widths
}

embedded_custom_widths_f = {}
embedded_custom_widths_00 = {"agents\n00\n"+k:custom_widths[k] for k in custom_widths if k not in ["diffusion", "multibody"]}
embedded_custom_widths_01 = {"agents\n01\n"+k:custom_widths[k] for k in custom_widths if k not in ["diffusion", "multibody"]}
embedded_custom_widths_f.update(embedded_custom_widths_00)
embedded_custom_widths_f.update(embedded_custom_widths_01)

config_grow_divide_topology_f = {
    'graph_format': 'hierarchy',
    'dashed_edges':True,
    'show_ports':False,
    'remove_nodes':[
        'agents\n00\nrates\ngrowth_rate', 'agents\n00\nrates\ngrowth_noise',
        'agents\n01\nrates\ngrowth_rate', 'agents\n01\nrates\ngrowth_noise',
        'dimensions'
    ],
    'process_colors':{
        'agents\n00\ndivide_condition':'tab:cyan',
        'agents\n00\ndivision':'tab:cyan',
        'agents\n00\nglobals_deriver':'tab:gray',
        'agents\n00\ngrowth':'tab:gray',
        
        'agents\n01\ndivide_condition':'tab:cyan',
        'agents\n01\ndivision':'tab:cyan',
        'agents\n01\nglobals_deriver':'tab:gray',
        'agents\n01\ngrowth':'tab:gray',
    },
    'store_colors':{
        'agents':'tab:cyan',
        'agents\n00':'tab:cyan',
        'agents\n00\nrates':'tab:gray',
        'agents\n00\nboundary':'tab:brown',
        'agents\n00\nboundary\nmass':'tab:brown',
        'agents\n00\nboundary\ndivide':'tab:brown',
        
        'agents\n01':'tab:cyan',
        'agents\n01\nrates':'tab:gray',
        'agents\n01\nboundary':'tab:brown',
        'agents\n01\nboundary\nmass':'tab:brown',
        'agents\n01\nboundary\ndivide':'tab:brown',
    },
    'coordinates':{
        
        'agents':(3*u,  u),
        
        'agents\n00\ndivision':(u, -3*u),
        'agents\n00\ndivide_condition':(2*u,  -3*u),
        'agents\n00':(0, 0),
        'agents\n00\nrates':(-1*u, -u),
        'agents\n00\nboundary':(0*u, -u),
        'agents\n00\nboundary\nmass':(1*u, -1.5*u),
        'agents\n00\nboundary\ndivide':(2*u, -1.5*u),
        'agents\n00\ngrowth':(-1*u, -3*u),
        'agents\n00\nglobals_deriver':(0*u, -3*u),
        
        'agents\n01\ndivision':(5.5*u, -3*u),
        'agents\n01\ndivide_condition':(6.5*u, -3*u),
        'agents\n01':(4.5*u, 0),
        'agents\n01\nrates':(3.5*u, -u),
        'agents\n01\nboundary':(4.5*u, -u),
        'agents\n01\nboundary\nmass':(5.5*u, -1.5*u),
        'agents\n01\nboundary\ndivide':(6.5*u, -1.5*u),
        'agents\n01\ngrowth':(3.5*u, -3*u),
        'agents\n01\nglobals_deriver':(4.5*u, -3*u)

    },
    'node_labels':{
        'agents\n00':"agent\n0",
        'agents\n00\nglobals_deriver':'agent 0\nglobals\nderiver',
        'agents\n00\ndivide_condition':'agent 0\ndivide\ncondition',
        'agents\n00\ndivision':'agent 0\ndivision',
        'agents\n00\nrates':'agent 0\nrates',
        'agents\n00\ngrowth':'agent 0\ngrowth',
        'agents\n00\ngloabls\nderiver':'agent 0\nglobals\nderiver',
        'agents\n00\nboundary':'agent 0\nboundary',
        'agents\n00\nboundary\nmass':'agent 0\nboundary\nmass',
        'agents\n00\nboundary\ndivide':'agent 0\nboundary\ndivide',
        
        'agents\n01':"agent\n1",
        'agents\n01\nglobals_deriver':'agent 1\nglobals\nderiver',
        'agents\n01\ndivide_condition':'agent 1\ndivide\ncondition',
        'agents\n01\ndivision':'agent 1\ndivision',
        'agents\n01\nrates':'agent 1\nrates',
        'agents\n01\ngrowth':'agent 1\ngrowth',
        'agents\n01\ngloabls\nderiver':'agent 1\nglobals\nderiver',
        'agents\n01\nboundary':'agent 1\nboundary',
        'agents\n01\nboundary\nmass':'agent 1\nboundary\nmass',
        'agents\n01\nboundary\ndivide':'agent 1\nboundary\ndivide',
    },
    'custom_widths':embedded_custom_widths_f
}

config_diffusion_topology = {
    'dashed_edges':True,
    'show_ports':False,
    'process_colors':{
        'diffusion_field':'tab:blue',
    },
    'store_colors':{
        'fields':'tab:blue',
        'agents':'tab:cyan',
        'dimensions':'tab:blue'
    },
    'node_labels':{
        'diffusion_field':'diffusion'
    },
    'custom_widths':{
        'diffusion_field':8
        }
}

config_lattice_topology = {
    'dashed_edges':True,
    'show_ports':False,
    'process_colors':{
        'multibody':'tab:blue',
        'diffusion':'tab:blue',
    },
    'store_colors':{
        'fields':'tab:blue',
        'agents':'tab:cyan',
        'dimensions':'tab:blue'
    },
    'coordinates': {
        'diffusion': (1.5, 1),
        'multibody': (2.5, 1),
    },
    'custom_widths':custom_widths
}

u = .8
config_grow_divide_lattice_topology = {
    'graph_format': 'hierarchy',
    'dashed_edges': True,
    'show_ports': False,
    'remove_nodes':[
        'agents\n0\nrates\ngrowth_rate', 'agents\n0\nrates\ngrowth_noise',
        # 'dimensions'
    ],
    'process_colors':{
        'agents\n0\ndivide_condition':'tab:cyan',
        'agents\n0\ndivision':'tab:cyan',
        'multibody':'tab:blue',
        'diffusion':'tab:blue',
        'agents\n0\nglobals_deriver':'tab:gray',
        'agents\n0\ngrowth':'tab:gray',
    },
    'store_colors':{
        'agents':'tab:cyan',
        'agents\n0':'tab:cyan',
        'agents\n0\nrates':'tab:gray',
        'agents\n0\nboundary':'tab:brown',
        'agents\n0\nboundary\nmass':'tab:brown',
        'agents\n0\nboundary\nexternal':'tab:brown',
        'agents\n0\nboundary\ndivide':'tab:brown',
    },
    'coordinates':{
        'agents\n0\ndivision':(1.5*u, -3*u),
        'agents\n0\ndivide_condition':(2.5*u, -3*u),
        'multibody':(2*u, 2*u),
        'diffusion':(3*u, 2*u),
        
        'agents':(0,  u),
        'agents\n0':(0, 0),
        'agents\n0\nrates':(-1*u, -u),
        'agents\n0\nboundary':(0*u, -u),
        'agents\n0\nboundary\nexternal':(1*u, -1.5*u),
        'agents\n0\nboundary\nmass':(2*u, -1.5*u),
        'agents\n0\nboundary\ndivide':(3*u, -1.5*u),
        
        'fields':(2.5*u, 1*u),
        'dimensions':(3.5*u, 1*u),
        'agents\n0\ngrowth':(-.5*u, -3*u),
        'agents\n0\nglobals_deriver':(.5*u, -3*u)
    },
    'node_labels':{
        'agents\n0':"agent\nn",
        'agents\n0\nglobals_deriver':'globals\nderiver',
        'agents\n0\ndivide_condition':'divide\ncondition',
        'agents\n0\ndivision':'division',
        'agents\n0\nrates':'rates',
        'agents\n0\ngrowth':'growth',
        'agents\n0\ngloabls\nderiver':'globals\nderiver',
        'agents\n0\nboundary':'boundary',
        'agents\n0\nboundary\nexternal':'boundary\nexternal',
        'agents\n0\nboundary\nmass':'boundary\nmass',
        'agents\n0\nboundary\ndivide':'boundary\ndivide',
    },
    'custom_widths':embedded_custom_widths
}


su = .8
RS = 0 #store row
RP = -1.8*su #process row

# node coordinates
agent_coordinates = {       
        # cobra
        'cobra':(-4*su, RP),
        'mass_deriver':(-3*su, RP),
        'volume_deriver':(-2*su, RP),
        'reactions':(-4*su, RS),
        'internal_counts':(-3*su, RS),
        'cobra_external':(-2*su, RS),
        'flux_bounds':(-1*su, RS),
        'fields':(4*su, RS-.5*su),

        # bioscrape
        'flux_adaptor':(-1*su, RP),
        'bioscrape':(0*su, RP),
        'dilution_rate_adaptor':(1*su, RP),
        'biomass_adaptor':(2*su, RP),
        'delta_counts_to_concs':(1*su, RP),
        'species':(0, RS),
        'rates':(su, RS),
        'delta_species':(2*su, RS),
        'delta_concs':(3*su, RS),
        
        # boundary/division
        'local_field':(5*su, RP),
        'field_counts_deriver':(6*su, RP),
        'division':(3*su, RP),
        'divide_condition':(4*su, RP),
        'boundary':(3*su, RS),
        'boundary\nexternal':(5*su, RS+.5*su),
        'boundary\nlocation':(4*su, RS+.5*su),
        'boundary\nexchanges':(5.5*su, RS-su),
        'boundary\nmass':(4*su, RS-.5*su),
        'boundary\ndivide':(5*su, RS-.5*su)}
global_coordinates = {
        'agents':(0, RS+2*su),
        'agents\n0':(0, RS+su),
        'multibody':(5.5*su, RS+2*su),
        'diffusion':(6.5*su, RS+2*su),
        'dimensions': (6*su, RS+su),
        'fields':(7*su, RS+su),
        # 'fields': (2.5 * su, 1 * su),
}
embedded_coordinates = {
    'agents\n0\n'+node_id: coord 
    for node_id, coord in agent_coordinates.items()}
embedded_coordinates.update(global_coordinates)



# Bioscrape Cobra Config
# node labels
agent_node_labels = {
        'flux_bounds':'flux\nbounds',
        'internal_counts':'internal\ncounts',
        'cobra_external':'cobra\nexternal',
        'hidden_exchanges':'hidden\nexchanges',
        'flux_adaptor':'flux\nadaptor',
        'dilution_rate_adaptor':'dilution rate\nadaptor',
        'biomass_adaptor':'biomass\nadaptor',
        'mass_deriver':'mass\nderiver',
        'volume_adaptor':'volume\nadaptor',
        'delta_species:':'delta\nspecies',
        'mass_deriver':'mass\nderiver',
        'volume_deriver':'volume\nderiver',
        'local_field': 'local\nfield',
        'boundary\nno_units':'boundary\n(no units)',
        'strip_units':'strip\nunits',
        'delta_concentrations':'delta\nconcs',
        'delta_counts_to_concs':'delta counts\nto\nconcs',
        'delta_species':'delta\nspecies',
        'field_counts_deriver':"field\ncounts\nderiver",

        'cobra':'cobra',
        'bioscrape':'bioscrape',
        'flux_adaptor':'flux\nadaptor',
        'dilution_rate_adaptor':'dilution rate\nadaptor',
        'clock':'clock',
        'mass_deriver':'mass\nderiver',
        'volume_deriver':'volume\nderiver',
        'biomass_adaptor':'biomass\nadaptor',
        'strip_units':'strip\nunits',
        'local_field':'local\nfield',
        'divide_condition':'divide\ncondition',
        'division':'division',
        'delta_counts_to_concs':'delta counts\nto\nconcs',
        'field_counts_deriver':'field\ncounts\nderiver',
        
        'species':'species',
        'internal_counts':"internal\ncounts",
        'delta_species':'delta\nspecies',
        'reactions':'reactions',
        'internal_counts':'internal\ncounts',
        'rates':'rates',
        'boundary':'boundary',
        'cobra_external':'cobra\nexternal',
        'hidden_exchanges':'hidden\nexchanges',
        'flux_bounds':'flux\nbounds',
        'boundary\ntime':'boundary\ntime',
        'boundary\nexternal':'boundary\nexternal',
        'boundary\ndivide':'boundary\ndivide',
        'boundary\nmass':'boundary\nmass',
        'boundary\nexchanges':'boundary\nexchanges',
        'boundary\nlocation':'boundary\nlocation',
        'boundary\nno_units':'boundary\n(no units)'
    }
global_node_labels = {
        'agents\n0':'agent\nn'}
embedded_node_labels = {
    'agents\n0\n'+node_id: label 
    for node_id, label in agent_node_labels.items()}
embedded_node_labels.update(global_node_labels)


# node colors
agent_process_colors = {
        'cobra':'tab:orange',
        'mass_deriver':'tab:orange',
        'volume_deriver':'tab:orange',
        'field_counts_deriver':'tab:brown',
        'biomass_adaptor':'tab:green',
        'bioscrape':'tab:green',
        'flux_adaptor':'tab:green',
        'dilution_rate_adaptor':'tab:green',
        'delta_counts_to_concs':'tab:green',
        'clock':'tab:gray',
        'strip_units':'tab:brown',
        'local_field':'tab:brown',
        'division':'tab:cyan',
        'divide_condition':'tab:cyan'}
agent_store_colors = {
        'flux_bounds':'tab:orange',
        'internal_counts':'tab:orange',
        'reactions':'tab:orange',
        'hidden_exchanges':'tab:orange',
        'cobra_external':'tab:orange',
        'rates':'tab:green',
        'species':'tab:green',
        'delta_species':'tab:green',
#         'delta\nconcs':'tab:green',
        'agents':'tab:cyan',
        'agents\n0':'tab:cyan',
        'boundary':'tab:brown',
        'boundary\nmass':'tab:brown',
        'boundary\ndivide':'tab:brown',
        'boundary\nexternal':'tab:brown',
        'boundary\nlocation':'tab:brown',
        'boundary\nexchanges':'tab:brown',
        'boundary\n(no units)':'tab:brown',
        'boundary\ntime':'tab:brown',
        'dimensions':'tab:brown',
#         'field_counts_deriver':'tab:brown',
        }
global_process_colors = {
        'diffusion':'tab:blue',
        'multibody':'tab:blue'}
global_store_colors = {
        'agents': 'tab:cyan',
        'agents\n0': 'tab:cyan',
        'fields':'tab:blue'}
embedded_process_colors = {
    'agents\n0\n'+node_id: color 
    for node_id, color in agent_process_colors.items()}
embedded_store_colors = {
    'agents\n0\n'+node_id: color 
    for node_id, color in agent_store_colors.items()}
embedded_process_colors.update(global_process_colors)
embedded_store_colors.update(global_store_colors)

embedded_custom_widths = {"agents\n0\n"+k:custom_widths[k] for k in custom_widths if k not in ["diffusion", "multibody"]}
embedded_custom_widths['diffusion'] = custom_widths['diffusion']
embedded_custom_widths['multibody'] = custom_widths['multibody']

config_bioscrape_cobra_topology = {
    'graph_format': 'hierarchy',
    'dashed_edges':True,
    'show_ports':False,
    'remove_nodes':[
        'boundary\ndimensions', 
        'boundary\nno_units',
        'boundary\ntime',
        'clock',
        'strip_units',
        'connect_external',
        'hidden_exchanges',
        'dimensions',
        'delta_concentrations',
        'boundary\nexchanges',
        
        'agents\n0\nboundary\ndimensions', 
        'agents\n0\nboundary\nno_units',
        'agents\n0\nboundary\ntime',
        'agents\n0\nclock',
        'agents\n0\nstrip_units',
        'agents\n0\nconnect_external',
        'agents\n0\nhidden_exchanges',
        'agents\n0\ndelta_concentrations',
        'agents\n0\nboundary\nexchanges'
    ],

}

config_single_cell_bioscrape_cobra_topology = copy.deepcopy(config_bioscrape_cobra_topology)
config_single_cell_bioscrape_cobra_topology['coordinates'] = agent_coordinates
config_single_cell_bioscrape_cobra_topology['node_labels'] = agent_node_labels
config_single_cell_bioscrape_cobra_topology['store_colors'] = agent_store_colors
config_single_cell_bioscrape_cobra_topology['process_colors'] = agent_process_colors
config_single_cell_bioscrape_cobra_topology['custom_widths'] = custom_widths

config_embedded_bioscrape_cobra_topology = copy.deepcopy(config_bioscrape_cobra_topology)
config_embedded_bioscrape_cobra_topology['coordinates'] = embedded_coordinates
config_embedded_bioscrape_cobra_topology['node_labels'] = embedded_node_labels
config_embedded_bioscrape_cobra_topology['process_colors'] = embedded_process_colors
config_embedded_bioscrape_cobra_topology['store_colors'] = embedded_store_colors
config_embedded_bioscrape_cobra_topology['custom_widths'] = embedded_custom_widths
config_embedded_bioscrape_cobra_topology['remove_nodes'].remove('dimensions')