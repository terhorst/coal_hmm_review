import numpy as np
import luigi
import collections
import sys

from config import *
import estimate.tasks

@luigi.util.inherits(estimate.tasks.PairwiseMomiAnalysisFromSimulatedData)
class BootstrapAnalysis(luigi.Task):
    def requires(self):
        return [
                self.clone(estimate.tasks.PairwiseMomiAnalysisFromSimulatedData)
                for _ in range(GlobalConfig().bootstrap_replicates)
                ]

    def run(self):
        for rep in self.input():
            result = unpickle(rep)
            print(result)
            aoeu

    def output(self):
        return GlobalConfig().local_target(
                "bootstrap", 
                self.seed, 
                "-".join(self.populations),
                "estimates.dat")

