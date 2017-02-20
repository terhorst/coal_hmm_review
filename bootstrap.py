import numpy as np
import luigi
import collections
import sys

from config import *
import data.tasks
import data.original
from estimate.tasks import PairwiseMomiAnalysis
from simulate.tasks import SimulateTwoPopulationSplitModel, MsprimeToVcf

class BootstrapSimulate(luigi.Task):
    '''
    Perform bootstrap analysis on fitted data set:

        1. Obtain MLE
        2. Simulate under MLE

    '''
    seed = luigi.IntParameter()
    populations = luigi.ListParameter()

    def requires(self):
        return data.tasks.PopulationMap()

    def run(self):
        np.random.seed(self.seed)
        population_map = unpickle(self.input())
        sample_names = [name for population in self.populations 
                             for name in population_map[population]]
        msprime_dumps = yield [
            SimulateTwoPopulationSplitModel(
                    seed=np.random.randint(2**32 - 1),
                    populations=self.populations)
            for _ in range(30)]
        vcfs = yield [MsprimeToVcf(
            msprime_dump_path=dump.path, 
            sample_names=sample_names)
            for dump in msprime_dumps]
        momi_datasets = yield [data.tasks.VCF2Momi(vcf_path=vcf.path) for vcf in vcfs]
        sfs = collections.Counter()
        for momi_dataset in momi_datasets:
            sfs.update(unpickle(momi_dataset)['sfs'])
        # Dump output
        self.output().makedirs()
        pickle.dump(sfs, open(self.output().path, "wb"), -1)

    def output(self):
        return GlobalConfig().local_target(
                "bootstrap", 
                self.seed, 
                "-".join(sorted(self.populations)),
                "simulations.dat")


class BootstrapEstimate(luigi.Task):
    seed = luigi.IntParameter()
    populations = luigi.ListParameter()

    def requires(self):
        return BootstrapSimulate(
                seed=self.seed, 
                populations=self.populations)

    def run(self):
        self.output().makedirs()
        estimates = yield PairwiseMomiAnalysis(
                populations=self.populations, 
                data_set_paths=[self.input().path])
        self.output().makedirs()
        pickle.dump(estimates, open(self.output().path, "wb"))

    def output(self):
        return GlobalConfig().local_target(
                "bootstrap", 
                self.seed, 
                "-".join(self.populations),
                "estimates.dat")

