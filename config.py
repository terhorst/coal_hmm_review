import importlib
import luigi
import pickle
import os

HPC = os.environ.get("HPC", False)
print("HPC mode:", HPC)
PSMC_PATH = os.environ.get("PSMC_PATH", "/tmp/coal_hmm_review/psmc")
MSMC_PATH = os.environ.get("MSMC_PATH", "/tmp/coal_hmm_review/msmc")

MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = MUTATION_RATE / 10.
N0 = 1e4
THETA = MUTATION_RATE * 2 * N0
GENERATION_TIME = 29

def unpickle(task):
    return pickle.load(open(task.path, "rb"))

class GlobalConfig(luigi.Config):
    chromosome_length = luigi.IntParameter()
    n_contigs = luigi.IntParameter(10)
    output_directory = luigi.Parameter()

    def local_target(self, *args):
        lt = luigi.LocalTarget(
            os.path.join(
                self.output_directory, 
                *[str(a) for a in args]))
        lt.makedirs()
        return lt
