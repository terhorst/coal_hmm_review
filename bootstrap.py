import collections
import itertools
import luigi
import numpy as np
import sys

from config import *
import tasks

@luigi.util.inherits(tasks.PairwiseMomiAnalysisFromSimulatedData)
class BootstrapAnalysis(luigi.Task):
    def requires(self):
        np.random.seed(self.seed)
        return [
                self.clone(
                    tasks.PairwiseMomiAnalysisFromSimulatedData,
                    seed=np.random.randint(2 ** 32))
                for _ in range(GlobalConfig().bootstrap_replicates)
                ]

    def run(self):
        for rep in self.input():
            result = unpickle(rep)
            print(result)

    def output(self):
        return GlobalConfig().local_target(
                "bootstrap", 
                self.seed, 
                "-".join(self.populations),
                "analysis.dat")


class BootstrapAllPairs(luigi.Task):
    seed = luigi.IntParameter()
    def requires(self):
        return tasks.PopulationMap()

    def run(self):
        np.random.seed(self.seed)
        pmap = unpickle(self.input())
        tasks = []
        for c in itertools.combinations(pmap, 2):
            populations = sorted(c)  
            tasks.append(BootstrapAnalysis(
                seed=np.random.randint(2**32 - 1), 
                populations=populations))
        yield tasks

    def complete(self):
        return False
