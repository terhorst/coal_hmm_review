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
        pprint(ev)
        if True:
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


class MsprimeToVcf(luigi.Task):
    msprime_simulation = luigi.TaskParameter()
    sample_names = luigi.Parameter()

    def _requires(self):
        yield GlobalConfig().local_target(
                self.input().path + ".vcf.gz.tbi")
        yield luigi.Task._requires(self)

    def requires(self):
        return self.msprime_simulation

    def output(self):
        return GlobalConfig().local_target(
                self.input().path + ".vcf.gz")

    def run(self):
        new_samples = ['msp_%d %s' % t for t in enumerate(self.sample_names)]
        tree_seq = msprime.load(self.input().path)
        vcf_path = self.output().path[:-3]
        with open(vcf_path, "wt") as f:  # omit the .gz
            tree_seq.write_vcf(f, 2)
        with contextlib.ExitStack() as stack:
            sample_renames = stack.enter_context(tempfile.NamedTemporaryFile("wt"))
            out_vcf_gz = stack.enter_context(open(self.output().path, 'wb'))
            open(sample_renames.name, "wt").write("\n".join(new_samples))
            bgzip(bcftools('reheader', '-s', 
                sample_renames.name, vcf_path, _piped=True), 
                    "-c", _out=out_vcf_gz)
        tabix(self.output().path)
