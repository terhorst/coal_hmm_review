import string
import msprime as msp
import attr
import numpy as np

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

    def to_csv(self, generation_time):
        return self.events.to_csv(generation_time)

    @classmethod
    def factory(cls, demo, N):
        events = DEMOGRAPHIES[demo]
        return cls(events, [N] * events.npop)

@attr.s
class DemographicEvents:
    npop = attr.ib()
    events = attr.ib()

    def __iter__(self):
        return iter(self.events)
    
    def to_csv(self, generation_time):
        out = []
        first = True
        for pid in range(self.npop):
            def pred(ev):
                if ev.type == "population_parameters_change":
                    return ev.population == pid
                elif ev.type == "mass_migration":
                    return ev.source == pid
                return False
            events_i = [ev for ev in self.events if pred(ev)]
            out.append(_events_to_csv(events_i, generation_time, first))
            first = False
        return "\n".join(out)
            

def _events_to_csv(events, generation_time, header=True):
    out = ["t,Ne,method"] * header
    ev = events[0]
    Ne = ev.initial_size
    g = ev.growth_rate
    t = ev.time
    tag = f"truth{ev.population}"
    out.append(f"{t * generation_time},{Ne},{tag}")
    for ev in events[1:]:
        if ev.type == 'population_parameters_change':
            if g != 0:
                tt = np.geomspace(max(t, 1), ev.time, 100, endpoint=False)
            else:
                tt = [ev.time]
            Ne_t = Ne * np.exp(-g * tt)
            out += [f"{ttt * generation_time},{Ne_tt},{tag}" for ttt, Ne_tt in zip(tt, Ne_t)]
            Ne = ev.initial_size
            g = ev.growth_rate
        t = ev.time
        out.append(f"{t * generation_time},{Ne},{tag}")
    if ev.type == "population_parameters_change":  # this assumes the joined-on lineage comes first
        out.append(f"{1e7},{Ne},{tag}")
    return "\n".join(out)

DEMOGRAPHIES = {
    'constant': DemographicEvents(1,
        [
            msp.PopulationParametersChange(
                time=0, initial_size=1e4, growth_rate=0, population_id=0),
            ]),

    'bottleneck': DemographicEvents(1, [
        msp.PopulationParametersChange(
            time=0, initial_size=1e5, growth_rate=0, population_id=0),
        msp.PopulationParametersChange(
            time=1e3, initial_size=5e3, population_id=0, growth_rate=0),
        msp.PopulationParametersChange(
            time=2e3, initial_size=1e5, population_id=0, growth_rate=0),
    ]),
    'recent_growth': DemographicEvents(1, [
        msp.PopulationParametersChange(
            time=0, initial_size=5e5, growth_rate=1e-2, population_id=0),
        msp.PopulationParametersChange(
            time=5e2, initial_size=2e4, population_id=0, growth_rate=0),
    ]),
    'migration': DemographicEvents(2, [
        msp.PopulationParametersChange(
            time=0, initial_size=1e6, growth_rate=0, population_id=0),
        msp.PopulationParametersChange(
            time=0, initial_size=2e4, population_id=1, growth_rate=0),
        msp.MigrationRateChange(time=500, rate=1/500, matrix_index=(0, 1)),
        msp.MigrationRateChange(time=500, rate=1/500, matrix_index=(1, 0)),
        msp.MassMigration(time=1000, source=0, dest=1, proportion=1.0)
    ])
}
