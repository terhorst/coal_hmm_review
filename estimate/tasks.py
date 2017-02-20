import os
import luigi
import pickle
import collections

import data.tasks
import estimate.momi
from config import *
from . import momi

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


class PairwiseMomiAnalysis(luigi.Task):
    populations = luigi.ListParameter()
    data_set_paths = luigi.ListParameter()

    def requires(self):
        return {'population_map': data.tasks.PopulationMap()}

    def build_sfs(self):
        sfs = collections.Counter()
        for path in self.data_set_paths:
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
