"""
Experiment to profile runtime in process next_update vs
remaining vivarium overhead
"""

import os
import cProfile
import pstats
from pstats import SortKey

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from vivarium.core.engine import Engine
from vivarium.core.control import run_library_cli
from vivarium.core.composition import EXPERIMENT_OUT_DIR

from bioscrape_cobra.bioscrape_cobra_stochastic import SBML_FILE_STOCHASTIC
from bioscrape_cobra.bioscrape_cobra_deterministic import SBML_FILE_DETERMINISTIC
from bioscrape_cobra.simulate import get_bioscrape_cobra_composite

# get the sbml files at the correct path
dirname = os.path.dirname(__file__)
DETERMINISTIC_FILE = os.path.join(dirname, SBML_FILE_DETERMINISTIC)
STOCHASTIC_FILE = os.path.join(dirname, SBML_FILE_STOCHASTIC)

DEFAULT_EXPERIMENT_TIME = 100
PROCESS_UPDATE_MARKER = 'bD'
VIVARIUM_OVERHEAD_MARKER = 'rD'
SIMULATION_TIME_MARKER = 'gD'


class ModelProfiler:
    """Profile Bioscrape-COBRA composites"""

    # model complexity
    n_agents = 1
    experiment_time = DEFAULT_EXPERIMENT_TIME
    parallel = False
    reuse_processes = False
    stochastic = False
    division = False
    spatial = False
    emit_step = 1

    # initialize
    composite = None
    experiment = None
    initial_state = None

    def set_parameters(
            self,
            n_agents=None,
            experiment_time=None,
            parallel=None,
            reuse_processes=None,
            emit_step=None,
            stochastic=None,
            division=None,
            spatial=None,
    ):
        self.n_agents = \
            n_agents if n_agents is not None else self.n_agents
        self.experiment_time = \
            experiment_time or self.experiment_time
        self.parallel = \
            parallel or self.parallel
        self.reuse_processes = \
            reuse_processes or self.reuse_processes
        self.emit_step = \
            emit_step or self.emit_step
        self.stochastic = \
            stochastic or self.stochastic
        self.division = \
            division or self.division
        self.spatial = \
            spatial or self.spatial

    def _generate_composite(self, **kwargs):
        initial_agent_states = [
            {'rates': {
                'k_leak': 0.005  # less leak -> less spontanteous expression
            }}
        ]

        self.composite, _, self.initial_state = get_bioscrape_cobra_composite(
            n_agents=self.n_agents,
            initial_agent_states=initial_agent_states,
            stochastic=self.stochastic,
            division=self.division,
            spatial=self.spatial,
            initial_glucose=1e1,
            initial_lactose=5e1,
            depth=0.5,
            diffusion_rate=2e-2,
            jitter_force=1e-5,
            bounds=[30, 30],
            n_bins=[30, 30],
            sbml_file=STOCHASTIC_FILE if self.stochastic else DETERMINISTIC_FILE,
            parallel=self.parallel,
            reuse_processes=self.reuse_processes,
        )

    def _initialize_experiment(self, **kwargs):
        self.experiment = Engine(
            processes=self.composite['processes'],
            topology=self.composite['topology'],
            initial_state=self.initial_state,
            **kwargs)

    def _run_experiment(self, **kwargs):
        self.experiment.update(kwargs['experiment_time'])
        self.experiment.end()

    def _get_emitter_data(self, **kwargs):
        _ = kwargs
        data = self.experiment.emitter.get_data()
        return data

    def _get_emitter_timeseries(self, **kwargs):
        _ = kwargs
        timeseries = self.experiment.emitter.get_timeseries()
        return timeseries

    def _profile_method(self, method, **kwargs):
        """The main profiling method and of the simulation steps

        Args
            method: the simulation step. For example self._run_experiment
        """
        profiler = cProfile.Profile()
        profiler.enable()
        method(**kwargs)
        profiler.disable()
        stats = pstats.Stats(profiler)
        return stats

    def run_profile(self):
        print('GENERATE COMPOSITE')
        self._profile_method(
            self._generate_composite)

        print('INITIALIZE EXPERIMENT')
        self._profile_method(
            self._initialize_experiment)

        print('RUN EXPERIMENT')
        self._profile_method(
            self._run_experiment, experiment_time=self.experiment_time)

        print('GET EMITTER DATA')
        self._profile_method(
            self._get_emitter_data)

    def profile_communication_latency(self):
        self._generate_composite()
        self._initialize_experiment(display_info=False)

        # profile the experiment
        stats = self._profile_method(
            self._run_experiment,
            experiment_time=self.experiment_time,
        )

        # get next_update runtime
        next_update_amount = ("next_update",)
        _, stats_list = stats.get_print_list(next_update_amount)

        process_update_time = 0
        for s in stats_list:
            process_update_time += stats.stats[s][3]

        # get total runtime
        experiment_time = stats.total_tt
        store_update_time = experiment_time - process_update_time

        # print_stats = stats.strip_dirs().sort_stats(-1).print_stats()
        # looping_stats = stats.sort_stats(SortKey.TIME).print_stats(20)

        return process_update_time, store_update_time


def run_scan(
        sim,
        scan_values=None,
):
    """Run a scan

    Args
        sim: the ModelProfiler object.
        scan_values: a list of dicts, with individual scan values.
    """
    scan_values = scan_values or []

    saved_stats = []
    for scan_dict in scan_values:
        # set the parameters
        sim.set_parameters(**scan_dict)

        print(f'{scan_dict}')

        # run experiment
        process_update_time, store_update_time = \
            sim.profile_communication_latency()

        # save data
        stat_dict = {
            **scan_dict,
            'process_update_time': process_update_time,
            'store_update_time': store_update_time,
        }
        saved_stats.append(stat_dict)

    return saved_stats


# Plotting functions
####################

def _make_axis(fig, grid, plot_n, patches, title='', label=''):
    ax = fig.add_subplot(grid[plot_n, 0])
    ax.set_xlabel(label)
    ax.set_ylabel('wall time (s)')
    ax.set_title(title)
    ax.legend(
        handles=patches,
        # ncol=2,
        # bbox_to_anchor=(0.5, 1.2),  # above
        bbox_to_anchor=(1.45, 0.65),  # to the right
    )
    return ax


def _get_patches(
        process=True,
        overhead=True,
        experiment=True,
):
    patches = []
    if process:
        patches.append(mpatches.Patch(
            color=PROCESS_UPDATE_MARKER[0],
            label="process updates"))
    if overhead:
        patches.append(mpatches.Patch(
            color=VIVARIUM_OVERHEAD_MARKER[0],
            label="vivarium overhead"))
    if experiment:
        patches.append(mpatches.Patch(
            color=SIMULATION_TIME_MARKER[0],
            label="simulation time"))
    return patches


def _add_stats_plot(
        ax,
        saved_stats,
        variable_name,
        process_update=False,
        vivarium_overhead=False,
        experiment_time=False,
        markersize=10,
):
    # plot saved states
    for stat in saved_stats:
        variable = stat[variable_name]
        process_update_time = stat['process_update_time']
        store_update_time = stat['store_update_time']

        if process_update:
            ax.plot(
                variable, process_update_time,
                PROCESS_UPDATE_MARKER, markersize=markersize)
        if vivarium_overhead:
            ax.plot(
                variable, store_update_time,
                VIVARIUM_OVERHEAD_MARKER, markersize=markersize)
        if experiment_time:
            experiment_time = process_update_time + store_update_time
            ax.plot(
                variable, experiment_time,
                SIMULATION_TIME_MARKER, markersize=markersize)


def plot_scan_results(
        saved_stats,
        plot_all=False,
        n_agents_plot=False,
        n_agents_update_only=False,
        parallel_plot=False,
        fig=None,
        grid=None,
        axis_number=None,
        row_height=3,
        title=None,
        out_dir=EXPERIMENT_OUT_DIR,
        filename='profile',
):
    axis_number = axis_number or 0
    plot_types = [
        n_agents_plot,
        n_agents_update_only,
        parallel_plot,
    ]

    if plot_all:
        n_agents_plot = True
        n_agents_update_only = True
        parallel_plot = True

    # make figure
    if fig:
        assert grid, "fig must provide grid for subplots"
    else:
        n_cols = 1
        n_rows = sum(plot_types)
        fig = plt.figure(figsize=(n_cols * 6, n_rows * row_height))
        grid = plt.GridSpec(n_rows, n_cols)

    # initialize axes
    if n_agents_plot:
        patches = _get_patches()
        ax = _make_axis(
            fig, grid, axis_number, patches, title,
            label='initial number of agents')

        _add_stats_plot(
            ax=ax, saved_stats=saved_stats,
            variable_name='n_agents',
            process_update=True,
            vivarium_overhead=True,
        )
        axis_number += 1

    if n_agents_update_only:
        patches = _get_patches()
        ax = _make_axis(
            fig, grid, axis_number, patches, title,
            label='initial number of agents')

        _add_stats_plot(
            ax=ax, saved_stats=saved_stats,
            variable_name='n_agents',
            process_update=True,
        )
        axis_number += 1

    if parallel_plot:
        patches = _get_patches()
        ax = _make_axis(
            fig, grid, axis_number, patches, title,
            label='parallel (False/True)')
        _add_stats_plot(
            ax=ax, saved_stats=saved_stats,
            variable_name='parallel',
            experiment_time=True)
        axis_number += 1

    # save
    if filename:
        plt.subplots_adjust(hspace=0.5)
        # plt.figtext(0, -0.1, filename, size=8)
        os.makedirs(out_dir, exist_ok=True)
        fig_path = os.path.join(out_dir, filename[0:100])
        fig.savefig(fig_path, bbox_inches='tight')
    return fig


# Individual scan functions
###########################

def scan_n_agents():
    n_agents = [n for n in range(1, 5)]
    scan_values = [{'n_agents': n} for n in n_agents]

    sim = ModelProfiler()
    sim.experiment_time = 60
    saved_stats = run_scan(sim,
                           scan_values=scan_values)
    plot_scan_results(saved_stats,
                      n_agents_plot=True,
                      n_agents_update_only=True,
                      filename=f'scan_n_agents_{max(n_agents)}')


def scan_parallel_processes():
    scan_values = [
        {'parallel': False},
        {'parallel': True},
    ]

    sim = ModelProfiler()
    sim.experiment_time = 20
    saved_stats = run_scan(sim,
                           scan_values=scan_values)
    plot_scan_results(saved_stats,
                      parallel_plot=True,
                      filename=f'scan_parallel')


def scan_agents_parallel():
    n_agents_scan = [n for n in range(1, 6, 1)]
    parallel_scan = [True, False]

    n_cols = 1
    n_rows = len(parallel_scan)
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 3))
    grid = plt.GridSpec(n_rows, n_cols)

    for axis_idx, parallel in enumerate(parallel_scan):
        scan_values = []
        for n_agents in n_agents_scan:
            scan_value = {
                'n_agents': n_agents,
                'parallel': parallel,
                'reuse_processes': False,
            }
            scan_values.append(scan_value)

        sim = ModelProfiler()
        saved_stats = run_scan(sim,
                               scan_values=scan_values)
        plot_scan_results(saved_stats,
                          n_agents_plot=True,
                          # parallel_plot=False,
                          row_height=2.5,
                          fig=fig,
                          grid=grid,
                          axis_number=axis_idx,
                          title=f'Parallel {parallel}',
                          )
    plot_scan_results({},
                      fig=fig,
                      grid=grid,
                      filename='scan_agents_parallel'
                      )


scans_library = {
    '0': scan_n_agents,
    '1': scan_parallel_processes,
    '2': scan_agents_parallel,
}

# python bioscrape_cobra/profile_runtime.py -n [name]
if __name__ == '__main__':
    run_library_cli(scans_library)
