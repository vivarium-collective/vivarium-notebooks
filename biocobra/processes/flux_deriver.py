from vivarium.core.process import Process


class FluxDeriver(Process):
    """ Bioscrape delta species to fluxes for constraining COBRA """

    defaults = {
        'time_step': 1.0,
        'flux_keys': {},  # key --> {option dictionary}
        'default_options': {  # default options if option_dictionary is empty
            "input_type": "delta",
            # "delta" corresponds to changes between timesteps. "amount" corresponds to the absolute quantity
        }
    }

    def __init__(self, config=None):
        super().__init__(config)

        self.prev_inputs = {k: None for k in self.parameters["flux_keys"]}

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

        input_type = self.parameters['flux_keys'][flux_key].get("input_type",
                                                                self.parameters['default_options']["input_type"])

        # Set delta
        if input_type == "delta":
            delta = inputs[flux_key]
        elif input_type == 'amount':
            if self.prev_inputs[flux_key] is None:
                delta = 0
            else:
                delta = inputs[flux_key] - self.prev_inputs[flux_key]
        else:
            raise ValueError(f"Unknown input_type: {input_type} for flux_key {flux_key}")

        return delta / dt

    def next_update(self, timestep, states):
        inputs = states["inputs"]

        update = {}
        update['fluxes'] = {}

        for flux_key in self.parameters['flux_keys']:
            # flux_dt = self.get_flux_interval(flux_key, timestep)
            # if flux_dt is not None:
            update['fluxes'][flux_key] = self.compute_flux(flux_key, timestep, inputs)

        self.prev_inputs = inputs
        return update


class DilutionFluxDeriver(FluxDeriver):
    # This is a non-negative "percentage" flux used to control a dilution rate of a CRN from a the rate of growth in biomass

    def next_update(self, timestep, states):
        inputs = states["inputs"]

        update = {}
        update['fluxes'] = {}

        for flux_key in self.parameters['flux_keys']:
            flux = self.compute_flux(flux_key, timestep, inputs)

            # Convert to a percent of the previous input
            if self.prev_inputs[flux_key] is None or (
                    self.prev_inputs[flux_key] == 0 and states["inputs"][flux_key] == 0):
                # Edge case: If both the current amount and the previous amount were 0, return 0
                update['fluxes'][flux_key] = 0

            elif self.prev_inputs[flux_key] == 0 and states["inputs"][flux_key] > 0:
                # Edge case: If both the previous amount was 0, us the current amount for the normalization
                update['fluxes'][flux_key] = flux / states["inputs"][flux_key]

            else:
                # Standard case
                update['fluxes'][flux_key] = flux / self.prev_inputs[flux_key]

            # Enforce non-negativity
            if update['fluxes'][flux_key] < 0:
                update['fluxes'][flux_key] = 0

        self.prev_inputs = inputs
        return update


class AverageFluxDeriver(FluxDeriver):
    # This is similar to a FluxDeriver, but fluxes are averaged between updates. This is useful to compute fluxes
    # which are used by a process with a larger dt than the processes controlling the flux, for example a stochastic
    # CRN may be updated frequently but FBA may be updated less frequently.
    pass
