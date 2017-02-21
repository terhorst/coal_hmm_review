import os
import pickle
import luigi

import sh
smc = sh.Command("smc++")
tabix = sh.Command("tabix")
bgzip = sh.Command("bgzip")
bcftools = sh.Command("bcftools")

MUTATION_RATE = RECOMBINATION_RATE = 1.25e-8
N0 = 1e4
THETA = MUTATION_RATE * 2 * N0

def unpickle(task):
    return pickle.load(open(task.path, "rb"))

class GlobalConfig(luigi.Config):
    contigs = list(map(str, range(1, 2)))
    bootstrap_replicates = luigi.IntParameter(20)
    chromosomes_per_bootstrap = luigi.IntParameter(30)
    chromosome_length = luigi.IntParameter(int(1e5))
    migration_density = luigi.FloatParameter(1e5 / 3.)
    input_directory = luigi.Parameter()
    output_directory = luigi.Parameter()

    def local_target(self, *args):
        return luigi.LocalTarget(
                os.path.join(
                    self.output_directory, *[str(a) for a in args]))
