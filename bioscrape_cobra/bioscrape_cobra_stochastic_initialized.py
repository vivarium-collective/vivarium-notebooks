from bioscrape_cobra.bioscrape_cobra_stochastic import (
    BioscrapeCOBRAstochastic, Bioscrape, DynamicFBA,
    TreeMass, Clock, MassToCount, CountsToMolar, MolarToCounts, StripUnits,
    DivideCondition, MetaDivision, Volume, LocalField, AverageFluxAdaptor)


# The stochastic Bioscrape/COBRA composer
class BioscrapeCOBRAstochastic_initialized(BioscrapeCOBRAstochastic):

    def __init__(self, config=None):
        super().__init__(config)
        self.processes_initialized = False

    def initialize_processes(self, config):
        # Processes
        self.cobra = DynamicFBA(config['cobra'])
        self.bioscrape = Bioscrape(config['bioscrape'])
        self.clock = Clock(config['clock'])
        # self.flux_adaptor = AverageFluxAdaptor(config['flux_adaptor'])

        # Derivers
        self.mass_deriver = TreeMass()
        self.volume_deriver = Volume()
        self.delta_counts_to_concs = CountsToMolar(config['delta_counts_to_concs'])
        self.biomass_adaptor = MassToCount(config['mass_to_counts'])
        self.strip_units = StripUnits(config['strip_units'])
        self.local_field = LocalField(config['local_fields'])
        self.field_counts_deriver = MolarToCounts(config['field_counts_deriver'])

        self.processes_initialized = True


    def generate_processes(self, config):
        if not self.processes_initialized:
            self.initialize_processes(config)

        processes = {
            # Processes
            'cobra': self.cobra,
            'bioscrape': self.bioscrape,
            'clock': self.clock,
            'flux_adaptor': AverageFluxAdaptor(config['flux_adaptor']),  # flux adaptor has internal state

            # Derivers
            'mass_deriver': self.mass_deriver,
            'volume_deriver': self.volume_deriver,
            'delta_counts_to_concs': self.delta_counts_to_concs,
            'biomass_adaptor': self.biomass_adaptor,
            'strip_units': self.strip_units,
            'local_field': self.local_field,
            'field_counts_deriver': self.field_counts_deriver}

        # Division Logic
        if config['divide_on']:
            # division config
            daughter_path = config['daughter_path']
            agent_id = config['agent_id']
            division_config = dict(
                config.get('division', {}),
                daughter_path=daughter_path,
                agent_id=agent_id,
                composer=self)
            processes.update({
                'divide_condition': DivideCondition(config['divide_condition']),
                'division': MetaDivision(division_config)})

        return processes
