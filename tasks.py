import collections
import contextlib
import hashlib
import itertools
import luigi
import luigi.util
import msprime
import numpy as np
import os
import pickle
import pysam
import re
import tempfile

from config import *
import data.original
import estimate.momi
import simulate.msprime
import util

# All luigi tasks have to go in one big file for now due to circular
# dependencies and luigi.util.require/inherit. Yay.

### DATA RELATED TASKS

## Derived data sets
class PopulationMap(luigi.Task):
    """
    Dict populations => samples. Restricted to have only samples / pops
    that are actually in the VCF.
    """
    def requires(self):
        return {"vcf": data.original.OriginalFullVCF(chromosome=1),
                "populations": data.original.OriginalPopulations()}

    def output(self):
        return GlobalConfig().local_target("vcf_population_map.dat")

    def run(self):
        pops = {}
        with pysam.VariantFile(self.input()['vcf'].path) as vcf:
            vcf_samples = set(vcf.header.samples)
        with open(self.input()['populations'].path, "rt") as f:
            fields = next(f)[1:].strip().split("\t")
            for line in f:
                # The provided spreadsheet is malformatted. We can
                # restrict to the first len(fields) entries and the
                # last will be the population ID, with Mbororo Fulani
                # truncated to Mbororo.
                record = dict(zip(fields,
                                  re.split(r"\t+", line.strip())[:len(fields)]
                                  ))
                inter = set(record.values()) & vcf_samples
                assert len(inter) < 2
                if len(inter) == 1:
                    pops.setdefault(record['Population'].replace(" ", "_"), []).append(next(iter(inter)))
                else:
                    print(record)
        pops = {pop: list(vcf_samples & set(samples)) for pop, samples in pops.items()}
        pops = {pop: samples for pop, samples in pops.items() if samples}
        pickle.dump(collections.OrderedDict(pops), open(self.output().path, "wb"), -1)


class _VCFConverter(luigi.Task):
    chromosome = luigi.IntParameter()

    def requires(self):
        return {"vcf": data.original.OriginalFullVCF(chromosome=self.chromosome),
                "centromeres": data.original.IndexedCentromeres(),
                "populations": PopulationMap()}
        return ret

    @property
    def populations(self):
        return pickle.load(open(self.input()['populations'].path, "rb"))


class VCF2SMC(_VCFConverter):
    population = luigi.Parameter()
    distinguished = luigi.Parameter()

    def output(self):
        return GlobalConfig().local_target(
                    "smc", "data",
                    self.population,
                    "{}.{}.txt.gz".format(self.distinguished, self.chromosome))

    def run(self):
        # Composite likelihood over first 3 individuals
        samples = self.populations[self.population]
        undistinguished = set(samples) - set([self.distinguished])
        self.output().makedirs()
        smc("vcf2smc",
                # "-m", self.input()['centromeres'].path,
                '-d', self.distinguished, self.distinguished,
                self.input()['vcf'].path,
                self.output().path,
                self.chromosome,
                "{}:{}".format(self.population, ",".join(samples)))


class ConvertAllSMC(luigi.Task):
    def completed(self):
        return False

    def requires(self):
        return PopulationMap()

    def run(self):
        pops = unpickle(self.input())
        tasks = []
        for i in range(1, 23):
            for pop in pops:
                for dist in pops[pop][:3]:
                    tasks.append(VCF2SMC(chromosome=i, population=pop, distinguished=dist))
        yield tasks

         
class _VCF2Momi(luigi.Task):
    contig = luigi.Parameter()
    populations = luigi.ListParameter()

    def requires(self):
        return {"population_map": PopulationMap()}

    def output(self):
        return GlobalConfig().local_target("vcf2momi", 
                self.contig, "-".join(self.populations), "momi.dat")

    def run(self):
        pmap = unpickle(self.input()['population_map'])
        sfs = collections.Counter()
        pops = list(self.populations)
        with pysam.VariantFile(self.input()['vcf'].path) as vcf:
            for record in vcf.fetch(contig=self.contig):
                d = {}
                for pop in pops:
                    gts = [x for sample in pmap[pop]
                           for x in record.samples[sample]['GT']
                           if x is not None]
                    n = len(gts)
                    a = sum(gts)
                    d[pop] = (n - a, a)
                k = tuple([d[pop] for pop in pops])
                sfs[k] += 1
        n = {pop: 2 * len(pmap[pop]) for pop in pops}
        self.output().makedirs()
        pickle.dump({'sfs': sfs, 'populations': pops, 'n': n}, 
                     open(self.output().path, "wb"), -1)

class _MergeVCFToMomi(luigi.Task):
    def run(self):
        sfs = collections.Counter()
        for result in self.input():
            d = unpickle(result)
            sfs.update(d['sfs'])
        d['sfs'] = sfs
        self.output().makedirs()
        pickle.dump(d, open(self.output().path, "wb"), -1)

@luigi.util.inherits(_VCF2Momi)
class _OriginalVCFToMomi(_VCF2Momi):
    def requires(self):
        ret = _VCF2Momi.requires(self)
        ret['vcf'] = self.clone(data.original.OriginalFullVCF)
        return ret

class OriginalVCFToMomi(_MergeVCFToMomi):
    populations = luigi.ListParameter()

    def requires(self):
        return [self.clone(_OriginalVCFToMomi, contig=str(i)) for i in range(1, 23)]

    def output(self):
        return GlobalConfig().local_target(
                "momi", 
                "original",
                "-".join(self.populations), 
                "momi.combined.dat")

### Estimation-related tasks
class EstimateAllSizeHistories(luigi.Task):
    def requires(self):
        return PopulationMap()

    def run(self):
        pops = pickle.load(open(self.input(), "rb"))
        yield [EstimateSizeHistory(population=pop) for pop in pops]

class EstimateSizeHistory(luigi.Task):
    population = luigi.Parameter()

    def requires(self):
        return PopulationMap()

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
                VCF2SMC(
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
        data_set = unpickle(self.input())
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
        return luigi.LocalTarget(
                self.input().path + '.analysis.dat')

    def run(self):
        sfs, n = self.build_sfs()
        n = [n[pop] for pop in self.populations]
        mle = estimate.momi.PairwiseMomiEstimator(self.populations, sfs, n)
        demography = mle.run()
        self.output().makedirs()
        pickle.dump(demography, open(self.output().path, "wb"), -1)

class PairwiseMomiAnalysisFromOriginalData(_PairwiseMomiAnalysis):
    def requires(self):
        return self.clone(OriginalVCFToMomi)

@luigi.util.requires(PairwiseMomiAnalysisFromOriginalData)
class SimulatePairFromMLE(luigi.Task):
    seed = luigi.IntParameter()

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

@luigi.util.inherits(SimulatePairFromMLE)
class MsprimeToVcf(luigi.Task):
    def requires(self):
        return {'population_map': PopulationMap(),
                'msprime.dat': self.clone(SimulatePairFromMLE)}

    def complete(self):
        if not os.path.exists(self.output().path + ".tbi"):
            return False
        return luigi.Task.complete(self)

    def output(self):
        return luigi.LocalTarget(self.input()['msprime.dat'].path + ".vcf.gz")

    @property
    def sample_names(self):
        pmap = unpickle(self.input()['population_map'])
        return [sample 
                for pop in self.populations 
                for sample in pmap[pop]]

    def run(self):
        new_samples = ['msp_%d %s' % t for t in enumerate(self.sample_names)]
        tree_seq = msprime.load(self.input()['msprime.dat'].path)
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

@luigi.util.inherits(MsprimeToVcf)
class SimulatedVCFToMomi(_VCF2Momi):
    def requires(self):
        ret = _VCF2Momi.requires(self)
        ret['vcf'] = self.clone(MsprimeToVcf)
        return ret

@luigi.util.inherits(SimulatedVCFToMomi)
class PairwiseMomiAnalysisFromSimulatedData(_PairwiseMomiAnalysis):
    def requires(self):
        np.random.seed(self.seed)
        return [self.clone(SimulatedVCFToMomi, seed=np.random.randint(2 ** 32)) 
                for _ in range(GlobalConfig().chromosomes_per_bootstrap)]


class _AllPairs(luigi.Task):
    def requires(self):
        return PopulationMap()

    @property
    def population_map(self):
        return unpickle(self.input())

    def all_pairs_tasks(self):
        pmap = unpickle(self.input())
        tasks = {}
        for c in itertools.combinations(pmap, 2):
            pair = tuple(sorted(c))
            tasks[pair] = PairwiseMomiAnalysisFromOriginalData(populations=pair)
        return tasks
