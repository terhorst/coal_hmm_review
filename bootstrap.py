import numpy as np
import luigi
import sys

from config import *
import data
from estimate.tasks import PairwiseMomiAnalysis
from simulate.tasks import SimulateTwoPopulationSplitModel
from simulate.msprime import MsprimeToVcf

class BootstrapSimulate(luigi.Task):
    '''
    Perform bootstrap analysis on fitted data set:

        1. Obtain MLE
        2. Simulate under MLE
        3. Refit

    '''
    seed = luigi.IntParameter()
    populations = luigi.ListParameter()

    def requires(self):
        return {'population_map': data.PopulationMap(),
                'momi_mle': PairwiseMomiAnalysis(
                    populations=self.populations,
                    data_provider=data.tasks.VCF2Momi(
                        vcf_provider=data.original.IndexedVCF()))}

    def run(self):
        np.random.seed(self.seed)
        population_map = pickle.load(open(self.input().path, "rb"))
        sample_names = [name for population in self.populations 
                             for name in population_map[population]]
        tasks = []
        for _ in range(30):
            t = SimulateTwoPopulationSplitModel(
                    seed=np.random.randint(2**32 - 1),
                    populations=self.populations)
            t = MsprimeToVcf(msprime_simulation=t, sample_names=sample_names)
            t = VCF2Momi(vcf_provider=t)
            tasks.append(t)
        momi_data = yield tasks
        sfs = collections.Counter()
        for f in momi_data:
            sfs.update(pickle.load(open(f.path, "rb")))
        pickle.dump(sfs, open(self.output().path, "wb"), -1)

class BootstrapEstimate(luigi.Task):
    seed = luigi.IntParameter()
    populations = luigi.ListParameter()

    def requires(self):
        return BootstrapSimulate(seed=seed, populations=populations)

    def run(self):
        
