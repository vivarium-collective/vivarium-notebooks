from vivarium.core.process import Deriver


class FluxDeriver(Deriver):
    """ Bioscrape delta species to fluxes for constraining COBRA """

    defaults = {
        'time_step': 1,
        'flux_keys': [],  # Used to compute fluxes from discrete changes (assumed to be deltas)
        'flux_options': {}  # key --> list of options: 'percent', 'amount', 'positive', 'ignore zeros'
    }

    def initial_state(self, config=None):
        return {}

    def ports_schema(self):
        return {
            'deltas': {
                flux_key: {
                    '_default': 0.0}
                for flux_key in self.parameters['flux_keys']
            },
            'amounts': {
                flux_key: {
                    '_default': 0.0,
                    '_updater': 'set'}
                for flux_key in self.parameters['flux_keys']
            },
            'fluxes': {
                flux_key: {
                    '_updater': 'set',
                    '_emit': True}
                for flux_key in self.parameters['flux_keys']
            }
        }

    def next_update(self, timestep, states):
        deltas = states['deltas']
        amounts = states['amounts']
        options = self.parameters['flux_options']
        dt = self.parameters['time_step']

        if not hasattr(self, 'prev_amounts'):
            self.prev_amounts = dict(amounts)

        update = {}
        update['fluxes'] = {}

        for key in self.parameters['flux_keys']:
            skip_update = False
            # Set delta
            if key not in options or 'amount' not in options[key]:
                delta = deltas[key]
            elif 'amount' in options[key]:
                delta = amounts[key] - self.prev_amounts[key]

            # Apply options
            # percent calculates flux as a percent of the previous amount
            if key in options and 'percent' in options[key]:
                if self.prev_amounts[key] == 0:  # to avoid dividing by 0
                    skip_update = True
                else:
                    delta = delta / (self.prev_amounts[key])

            # positive sets all negative fluxes to 0
            if key in options and 'positive' in options[key] and delta < 0:
                delta = 0

            # if delta is 0, does not update this flux
            # this is a hack which allows for variable update times between Cobra and Bioscrape
            if key in options and 'ignore zeros' in options[key] and delta == 0:
                skip_update = True

            if not skip_update:
                update['fluxes'][key] = delta / dt

        self.prev_amounts = amounts
        return update