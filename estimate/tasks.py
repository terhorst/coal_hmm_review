import os
import luigi
import pickle
import collections

import data
from config import *
from . import momi

class EstimateAllSizeHistories(luigi.Task):
    def requires(self):
        return data.PopulationMap()

    def run(self):
        pops = pickle.load(open(self.input(), "rb"))
        yield [EstimateSizeHistory(population=pop) for pop in pops]

class EstimateSizeHistory(luigi.Task):
    population = luigi.Parameter()

    def requires(self):
        return data.PopulationMap()

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
                data.VCF2SMC(
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


class PairwiseMomiAnalysis(luigi.Task):
    populations = luigi.ListParameter()
    data_provider = luigi.TaskParameter()

    def requires(self):
        return {'population_map': data.PopulationMap(),
                'sfs': self.data_provider}

    def build_sfs(self):
        sfs = collections.Counter()
        n = None
        for f in self.input()['sfss']:
            d = pickle.load(open(f.path, "rb"))
            if n is None:
                n = d['n']
            else:
                assert n == d['n']
            i = [d['populations'].index(p) for p in self.populations]
            sfs[None] = d['sfs'][None]
            del d['sfs'][None]
            for entry in d['sfs']:
                key = (entry[i[0]], entry[i[1]])
                if ((key[0][0] == key[1][0] == 0) or 
                    (key[0][1] == key[1][1] == 0)):
                    # Recode sites which are monomorphic in the subsample
                    key = None
                sfs[key] += d['sfs'][entry]
        return sfs, n

    def output(self):
        return GlobalConfig().local_target(
                "momi", "estimates", "-".join(self.populations) + ".dat")

    def run(self):
        self.output().makedirs()
        sfs, n = self.build_sfs()
        n = [n[pop] for pop in self.populations]
        mle = momi.PairwiseMomiEstimator(self.populations, sfs, n)
        res = mle.run()
        mle_events = mle.events(*res.x)
        pickle.dump(
                {'populations': list(self.populations), 'mle_events': mle_events, 'n': n},
                open(self.output().path, "wb"),
                -1)
