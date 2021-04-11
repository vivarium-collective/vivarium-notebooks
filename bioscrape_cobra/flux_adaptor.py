from vivarium.core.process import Process
from vivarium.library.units import Quantity


class FluxAdaptor(Process):
    """ Bioscrape delta species to fluxes for constraining COBRA """

    defaults = {
        'time_step': 1.0,
        'flux_keys': {},  # key --> {option dictionary}
        'default_options': {  # default options if option_dictionary is empty
            'input_type': 'delta',
            # 'delta' corresponds to changes between timesteps. 'amount' corresponds to the absolute quantity
        }
    }

    def __init__(self, config=None):
        super().__init__(config)

        self.prev_inputs = {k: None for k in self.parameters['flux_keys']}

    def initial_state(self, config=None):
        return {}

    def ports_schema(self):
        return {
            'inputs': {
                flux_key: {
                    '_default': 0.0}
                for flux_key in self.parameters['flux_keys']
            },
            'fluxes': {
                flux_key: {
                    '_updater': 'set',
                    '_emit': True}
                for flux_key in self.parameters['flux_keys']
            }
        }

    def compute_flux(self, flux_key, dt, inputs):
        # computes the flux for a specific flux key of the time interval dt

        input_type = self.parameters['flux_keys'][flux_key].get(
            'input_type', self.parameters['default_options']['input_type'])

        # Set delta
        if input_type == 'delta':
            delta = inputs[flux_key]
        elif input_type == 'amount':
            if self.prev_inputs[flux_key] is None:
                delta = 0
            else:
                delta = inputs[flux_key] - self.prev_inputs[flux_key]
        else:
            raise ValueError(f'Unknown input_type: {input_type} for flux_key {flux_key}')

        return delta / dt

    def next_update(self, timestep, states):
        inputs = states['inputs']

        update = {}
        update['fluxes'] = {}

        for flux_key in self.parameters['flux_keys']:
            # flux_dt = self.get_flux_interval(flux_key, timestep)
            # if flux_dt is not None:
            update['fluxes'][flux_key] = self.compute_flux(flux_key, timestep, inputs)

        self.prev_inputs = inputs
        return update


class DilutionFluxAdaptor(FluxAdaptor):
    """
    This is a non-negative "percentage" flux used to control a dilution rate of a CRN from a the rate of growth in biomass
    """
    def next_update(self, timestep, states):
        inputs = states['inputs']

        update = {}
        update['fluxes'] = {}

        for flux_key in self.parameters['flux_keys']:
            flux = self.compute_flux(flux_key, timestep, inputs)

            # Convert to a percent of the previous input
            if self.prev_inputs[flux_key] is None or (
                    self.prev_inputs[flux_key] == 0 and states['inputs'][flux_key] == 0):
                # Edge case: If both the current amount and the previous amount were 0, don't update the value
                pass
                #update['fluxes'][flux_key] = 0

            elif self.prev_inputs[flux_key] == 0 and states['inputs'][flux_key] > 0:
                # Edge case: If both the previous amount was 0, don't update the value
                pass
                #update['fluxes'][flux_key] = flux / states['inputs'][flux_key]

            else:
                # Standard case
                percent = flux / self.prev_inputs[flux_key]
                if isinstance(percent, Quantity):
                    percent = percent.magnitude
                update['fluxes'][flux_key] = percent

            # Enforce non-negativity
            if flux_key in update['fluxes'] and update['fluxes'][flux_key] < 0:
                update['fluxes'][flux_key] = 0

        self.prev_inputs = inputs
        return update


class AverageFluxAdaptor(FluxAdaptor):
    """
    This is similar to a FluxAdaptor, but fluxes are averaged between updates. This is useful to compute fluxes
    which are used by a process with a larger dt than the processes controlling the flux, for example a stochastic
    CRN may be updated frequently but FBA may be updated less frequently.
    """

    defaults = {
        'time_step': 1.0,
        'flux_keys': {},  # key --> {option dictionary}
        'default_options': {  # default options if option_dictionary is empty
            'input_type': 'delta',
            # 'delta' corresponds to changes between timesteps. 'amount' corresponds to the absolute quantity
            'window_size': 10  # Corresponds to the number of timesteps to average over
        }
    }

    def __init__(self, config=None):
        super().__init__(config)

        # Stores a list of previously computed fluxes
        self.prev_fluxes = {k: [] for k in self.parameters['flux_keys']}

    def next_update(self, timestep, states):
        inputs = states['inputs']

        update = {}
        update['fluxes'] = {}

        for flux_key in self.parameters['flux_keys']:
            window_size = self.parameters['flux_keys'][flux_key].get(
                'window_size', self.parameters['default_options']['window_size'])
            # flux_dt = self.get_flux_interval(flux_key, timestep)
            # if flux_dt is not None:
            flux = self.compute_flux(flux_key, timestep, inputs)
            self.prev_fluxes[flux_key].append(flux)

            # Compute the average fluxes
            if len(self.prev_fluxes[flux_key]) < window_size:
                update['fluxes'][flux_key] = sum(self.prev_fluxes[flux_key]) / len(self.prev_fluxes[flux_key])
            else:
                update['fluxes'][flux_key] = sum(self.prev_fluxes[flux_key][-window_size:]) / window_size

        self.prev_inputs = inputs

        return update
