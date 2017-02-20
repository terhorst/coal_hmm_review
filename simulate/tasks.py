import pickle
import luigi

from config import *
import estimate.tasks
from . import msprime

class SimulateTwoPopulationSplitModel(luigi.Task):
    seed = luigi.IntParameter()
    populations = luigi.ListParameter()

    def requires(self):
        return estimate.tasks.PairwiseMomiAnalysis(self.populations)

    def output(self):
        return GlobalConfig().local_target("simulation", 
                "-".join(sorted(self.populations)),
                "msprime." + str(self.seed) + ".dat")

    def run(self):
        self.output().makedirs()
        d = pickle.load(open(self.input().path, "rb"))
        sim = msprime.MsprimeMomiSimulator(d['populations'], d['n'], d['mle_events'])
        sim.run(self.seed, self.output().path)
