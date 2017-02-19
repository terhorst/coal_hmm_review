import os
import luigi

import convert_data
from config import *

class EstimateSizeHistory(luigi.Task):
    population = luigi.Parameter()

    def requires(self):
        return [convert_data.VCF2SMC(contig=str(c), population=self.population) for c in range(1, 23)]

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
        smc('estimate', '--theta', .00025, '-v', '-o', self._output_directory)

class PairwiseMomiAnalysis(luigi.Task):
    populations = luigi.ListParameter()

    def requires(self):
        return [convert_data.VCF2Momi(contig=c) for c in range(1, 23)]

    def build_splits(self):
        c = collections.Counter()
        for f in self.input():
            d = pickle.load(open(f.path, "rb"))
            i = [d['pops'].index(p) for p in self.pops]
            c[None] = d['sfs'][None]
            del d['sfs'][None]
            for entry in d['sfs']:
                key = (entry[i[0]], entry[i[1]])
                if ((key[0][0] == key[1][0] == 0) or 
                    (key[0][1] == key[1][1] == 0)):
                    # Recode sites which are monomorphic in the subsample
                    key = None
                c[key] += d['sfs'][entry]

    def run(self):
        pass
