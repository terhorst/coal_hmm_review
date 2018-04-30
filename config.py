import os
import pickle
import luigi
import importlib

MUTATION_RATE = RECOMBINATION_RATE = 1.25e-8
N0 = 1e4
THETA = MUTATION_RATE * 2 * N0

def unpickle(task):
    return pickle.load(open(task.path, "rb"))

class GlobalConfig(luigi.Config):
    chromosome_length = luigi.IntParameter(int(1e7))
    n_contigs = luigi.IntParameter(10)
    output_directory = luigi.Parameter()

    def local_target(self, *args):
        lt = luigi.LocalTarget(
            os.path.join(
                self.output_directory, 
                *[str(a) for a in args]))
        lt.makedirs()
        return lt
