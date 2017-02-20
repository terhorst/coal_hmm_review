from pprint import pprint
import msprime
import contextlib
import tempfile

from config import *

class MsprimeMomiSimulator:
    def __init__(self, populations, n, momi_events):
        self.populations = populations
        self.momi_events = momi_events
        self.n = n

    def events(self):
        msp_events = []
        for ty, t, *ev in self.momi_events:
            t *= 2. * N0  # msprime scale is, sensibly, generations
            if ty in ("-ep", "-ej"):
                if ty == "-ej":
                    p = 1.0
                else:
                    p = ev[2]
                msp_events.append(msprime.MassMigration(time=t, 
                    source=self.populations.index(ev[0]),
                    destination=self.populations.index(ev[1]), 
                    proportion=p))
            elif ty == "-en":
                msp_events.append(msprime.PopulationParametersChange(time=t, 
                    population_id=self.populations.index(ev[0]),
                    initial_size=2 * N0 * ev[1], growth_rate=0.))
        msp_events = sorted(msp_events, key=lambda ev: ev.time)
        return msp_events

    def run(self, seed, output_path, length=int(1e5)):
        pc = [msprime.PopulationConfiguration(sample_size=n) for n in self.n]
        ev = self.events()
        if False:
            msprime.DemographyDebugger(Ne=N0,
                    population_configurations=pc,
                    demographic_events=ev).print_history()
        sim = msprime.simulate(
                mutation_rate=MUTATION_RATE,
                recombination_rate=RECOMBINATION_RATE,
                demographic_events=ev,
                population_configurations=pc,
                random_seed=seed,
                Ne=N0,
                length=length)
        sim.dump(output_path)
