import numpy as np
import luigi
import collections
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

