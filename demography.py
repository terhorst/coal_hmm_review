import string
import msprime as msp
import attr

@attr.s
class Demography:
    events = attr.ib()
    N = attr.ib()

    def _sample_name_iter(self):
        i = s = 0
        while i < len(self.N):
            yield ['msp_%d' % (s + j) for j in range(self.N[i] // 2)]
            s += self.N[i] // 2
            i += 1

    def samples(self):
        return list(self._sample_name_iter())

    def populations(self):
        return {pop: sample_names for pop, sample_names in zip(
            string.ascii_lowercase, self._sample_name_iter(self.N))}

    def population_configs(self):
        return [msp.PopulationConfiguration(sample_size=nn) for nn in self.N]


DEMOGRAPHIES = {
    'constant': [],
    'bottleneck': [
        msp.PopulationParametersChange(
            time=0, initial_size=1e5, growth_rate=0, population_id=0),
        msp.PopulationParametersChange(
            time=1e3, initial_size=5e3, population_id=0, growth_rate=0),
        msp.PopulationParametersChange(
            time=2e3, initial_size=1e5, population_id=0, growth_rate=0),
    ],
    'recent_growth': [
        msp.PopulationParametersChange(
            time=0, initial_size=5e5, growth_rate=1e-2, population_id=0),
        msp.PopulationParametersChange(
            time=5e2, initial_size=2e4, population_id=0, growth_rate=0),
    ]
}

def factory(demo, N):
    return Demography(DEMOGRAPHIES[demo], [N])
