import collections
import contextlib
import luigi
import luigi.util
import os
import pickle
import tempfile

from config import *
import data.tasks
import estimate.momi
import simulate.msprime

class EstimateAllSizeHistories(luigi.Task):
    def requires(self):
        return data.tasks.PopulationMap()

    def run(self):
        pops = pickle.load(open(self.input(), "rb"))
        yield [EstimateSizeHistory(population=pop) for pop in pops]

class EstimateSizeHistory(luigi.Task):
    population = luigi.Parameter()

    def requires(self):
        return data.tasks.PopulationMap()

    @property
    def population_map(self):
        return pickle.load(open(self.input().path, "rb"))

    @property
    def _output_directory(self):
        return os.path.join(
                GlobalConfig().output_directory,
                "smc",
                "estimates",
                self.population)

    def output(self):
        return luigi.LocalTarget(
                os.path.join(self._output_directory, "model.final.json")
                )

    def run(self):
        # Create data sets from composite likelihood
        samples = self.population_map[self.population]
        smc_data_files = yield [
                data.tasks.VCF2SMC(
                    contig=str(c), 
                    population=self.population, 
                    distinguished=s) 
                for c in GlobalConfig().contigs
                for s in list(samples)[:3]]
        smc('estimate', 
                '--theta', .00025, '--reg', 10, '--blocks', 1, 
                "--knots", 20,
                '--no-initialize',
                '-v', '-o', self._output_directory,
                *[f.path for f in smc_data_files])


class _PairwiseMomiAnalysis(luigi.Task):
    populations = luigi.ListParameter()

    def build_sfs(self):
        sfs = collections.Counter()
        for path in self.momi_datasets:
            data_set = pickle.load(open(path, "rb"))
            n = data_set['n']
            i = [data_set['populations'].index(p) for p in self.populations]
            for entry in data_set['sfs']:
                key = (entry[i[0]], entry[i[1]])
                if ((key[0][0] == key[1][0] == 0) or 
                    (key[0][1] == key[1][1] == 0)):
                    continue
                sfs[key] += data_set['sfs'][entry]
        return sfs, n

    def output(self):
        return GlobalConfig().local_target(
                "momi", "estimates", "-".join(self.populations) + ".dat")

    def run(self):
        self.output().makedirs()
        sfs, n = self.build_sfs()
        n = [n[pop] for pop in self.populations]
        mle = estimate.momi.PairwiseMomiEstimator(self.populations, sfs, n)
        res = mle.run()
        mle_events = mle.events(*res.x)
        pickle.dump(
                {'populations': list(self.populations), 
                 'mle_events': mle_events, 'n': n},
                open(self.output().path, "wb"),
                -1)

class PairwiseMomiAnalysisFromOriginalData(_PairwiseMomiAnalysis):
    def requires(self):
        return [data.tasks.OriginalVCF2Momi()]

@luigi.util.inherits(data.tasks.SimulatedVCFToMomi)
class PairwiseMomiAnalysisFromSimulatedData(_PairwiseMomiAnalysis):
    def requires(self):
        return [self.clone(data.tasks.SimulatedVCFToMomi) 
                for _ in range(GlobalConfig().chromosomes_per_boostrap)]

# Tasks that should be grouped under simulate/, but have an estimation
# component as well. Basically they are here to get out of circular
# dependency hell.
@luigi.util.requires(PairwiseMomiAnalysisFromOriginalData)
class SimulatePairFromMLE(luigi.Task):
    def output(self):
        return GlobalConfig().local_target(
                "simulation", 
                "-".join(sorted(self.populations)),
                "msprime." + str(self.seed) + ".dat")

    def run(self):
        mle = unpickle(self.input())
        self.output().makedirs()
        sim = simulate.msprime.MsprimeMomiSimulator(
                self.populations, mle.n, mle.events)
        sim.run(self.seed, self.output().path)

@luigi.util.requires(SimulatePairFromMLE)
class MsprimeToVcf(luigi.Task):
    def requires(self):
        return data.tasks.PopulationMap()

    def complete(self):
        if not os.path.exists(self.output().path + ".tbi"):
            return False
        return luigi.Task.complete(self)

    def output(self):
        return luigi.LocalTarget(self.input().path + ".vcf.gz")

    @property
    def sample_names(self):
        pmap = unpickle(self.input())
        return [sample 
                for pop in self.populations 
                for sample in pmap[pop]]

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
