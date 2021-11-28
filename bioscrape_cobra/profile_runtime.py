"""
Experiment to profile runtime in process next_update vs
remaining vivarium overhead
"""

import os
import cProfile
import pstats

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from vivarium.core.engine import Engine
from vivarium.core.control import run_library_cli

from bioscrape_cobra.bioscrape_cobra_stochastic import SBML_FILE_STOCHASTIC
from bioscrape_cobra.bioscrape_cobra_deterministic import SBML_FILE_DETERMINISTIC
from bioscrape_cobra.simulate import get_bioscrape_cobra_composite


# get the sbml files at the correct path
dirname = os.path.dirname(__file__)
DETERMINISTIC_FILE = os.path.join(dirname, SBML_FILE_DETERMINISTIC)
STOCHASTIC_FILE = os.path.join(dirname, SBML_FILE_STOCHASTIC)

DEFAULT_EXPERIMENT_TIME = 100
PROCESS_UPDATE_MARKER = 'b.'
VIVARIUM_OVERHEAD_MARKER = 'r.'
SIMULATION_TIME_MARKER = 'g.'


class ModelProfiler:
    """Profile Bioscrape-COBRA composites"""

    # model complexity
    experiment_time = DEFAULT_EXPERIMENT_TIME

    # initialize
    composite = None
    experiment = None

    def set_parameters(
            self,
            experiment_time=None,
    ):
        self.experiment_time = \
            experiment_time or self.experiment_time

    def _generate_composite(self, **kwargs):

        parallel = True
        bounds = [30, 30]
        n_bins = [30, 30]
        depth = 0.5
        initial_agent_states = [
            {'rates': {
                'k_leak': 0.005  # less leak -> less spontanteous expression
            }}
        ]

        self.composite, _, _ = get_bioscrape_cobra_composite(
            initial_agent_states=initial_agent_states,
            stochastic=True,
            division=True,
            spatial=True,
            initial_glucose=1e1,
            initial_lactose=5e1,
            depth=depth,
            diffusion_rate=2e-2,
            jitter_force=1e-5,
            bounds=bounds,
            n_bins=n_bins,
            sbml_file=STOCHASTIC_FILE,
            parallel=parallel
        )

    def _initialize_experiment(self, **kwargs):
        self.experiment = Engine(
            processes=self.composite['processes'],
            topology=self.composite['topology'],
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
        cc, nc, tt, ct, callers = stats.stats[stats_list[0]]
        _ = cc
        _ = nc
        _ = tt
        _ = callers
        process_update_time = ct

        # get total runtime
        experiment_time = stats.total_tt

        # analyze
        store_update_time = experiment_time - process_update_time

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
        n_parallel_processes = scan_dict.get('number_of_parallel_processes', 0)

        # set the parameters
        sim.set_parameters(
            # number_of_parallel_processes=n_parallel_processes,
        )

        print(
            f'number_of_parallel_processes={n_parallel_processes} '
        )

        # run experiment
        process_update_time, store_update_time = sim.profile_communication_latency()

        # save data
        stat_dict = {
            'number_of_parallel_processes': n_parallel_processes,
            'process_update_time': process_update_time,
            'store_update_time': store_update_time,
        }
        saved_stats.append(stat_dict)

    return saved_stats


def make_axis(fig, grid, plot_n, patches, label=''):
    ax = fig.add_subplot(grid[plot_n, 0])
    ax.set_xlabel(label)
    ax.set_ylabel('runtime (s)')
    ax.legend(
        loc='upper center',
        handles=patches, ncol=2,
        bbox_to_anchor=(0.5, 1.2), )
    return ax


def get_patches(
        process=True,
        overhead=True,
        experiment=False
):
    patches = []
    if process:
        patches.append(mpatches.Patch(
            color=PROCESS_UPDATE_MARKER[0], label="process updates"))
    if overhead:
        patches.append(mpatches.Patch(
            color=VIVARIUM_OVERHEAD_MARKER[0], label="vivarium overhead"))
    if experiment:
        patches.append(mpatches.Patch(
            color=SIMULATION_TIME_MARKER[0], label="simulation time"))
    return patches


def plot_scan_results(
        saved_stats,
        plot_all=True,
        parallel_plot=False,
        out_dir='out/experiments',
        filename='profile',
):
    if plot_all:
        parallel_plot = True

    n_cols = 1
    n_rows = sum([
        parallel_plot,
    ])

    # make figure
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 3))
    grid = plt.GridSpec(n_rows, n_cols)

    # initialize axes
    plot_n = 0
    if parallel_plot:
        patches = get_patches(
            process=False,
            overhead=False,
            experiment=True)
        ax_depth = make_axis(
            fig, grid, plot_n, patches,
            label='number of parallel processes')
        plot_n += 1

    # plot saved states
    for stat in saved_stats:
        n_parallel_processes = stat['number_of_parallel_processes']
        process_update_time = stat['process_update_time']
        store_update_time = stat['store_update_time']

        if parallel_plot:
            experiment_time = process_update_time + store_update_time
            ax_depth.plot(
                n_parallel_processes, experiment_time, SIMULATION_TIME_MARKER)

    # adjustments
    plt.subplots_adjust(hspace=0.5)
    plt.figtext(0, -0.1, filename, size=8)

    # save
    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, filename[0:100])
    fig.savefig(fig_path, bbox_inches='tight')


# scan functions
def scan_parallel_processes():
    total_processes = 20
    n_parallel_processes = [i*2 for i in range(int(total_processes/2))]
    scan_values = [
        {
            'number_of_processes': total_processes,
            'number_of_parallel_processes': n
        } for n in n_parallel_processes
    ]

    sim = ModelProfiler()
    sim.experiment_time = 100
    saved_stats = run_scan(sim,
                           scan_values=scan_values)
    plot_scan_results(saved_stats,
                      plot_all=False,
                      parallel_plot=True,
                      filename=f'scan_parallel_processes_{total_processes}')


scans_library = {
    '0': scan_parallel_processes,
}

# python bioscrape_cobra/profile_runtime.py -n [name]
if __name__ == '__main__':
    run_library_cli(scans_library)
