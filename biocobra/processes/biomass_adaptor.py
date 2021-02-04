from vivarium.core.process import Deriver
from vivarium.library.units import units


class BiomassAdaptor(Deriver):
    """ Adapts COBRA mass variable to Bioscrape biomass """

    defaults = {
        'input_mass_units': 1.0 * units.fg,
        'input_volume_units': 1.0 * units.fL,
        'output_mass_units': 1.0 * units.mmolar,
        'output_volume_units': 1.0 * units.fL,
        'mass_species_molecular_weight': 1.0 * units.fg / units.molec
    }

    def initial_state(self, config=None):
        return {}

    def ports_schema(self):
        return {
            'input': {
                'mass': {
                    '_default': 1.0 * units.fg,
                },
                'volume': {
                    '_default': 1.0 * units.fL}
            },
            'output': {
                # the value used by Bioscrape
                'biomass': {
                    '_default': 1.0,
                    '_update': 'set',
                }
            }
        }

    def next_update(self, timestep, states):
        mass = states['input']['mass']

        # do conversion
        # Concentration = mass/molecular_weight/characteristic volume
        # Note: Biomass is also used to set Volume, so here we just set the scale
        mass_species_conc = mass / self.config['mass_species_molecular_weight'] / (
                    1 * self.config['output_volume_units'])

        update = {
            'output': {
                # Return the correct units, with units stripped away for bioscrape
                'biomass': mass_species_conc.to(self.config['output_mass_units']).magnitude
            }
        }
        return update
