import os
import luigi

import sh
smc = sh.Command("smc++")
tabix = sh.Command("tabix")
bgzip = sh.Command("bgzip")

MUTATION_RATE = RECOMBINATION_RATE = 1.25e-8
N0 = 1e4
THETA = MUTATION_RATE * 2 * N0

class GlobalConfig(luigi.Config):
    contigs = list(map(str, range(1, 2)))
    input_directory = luigi.Parameter()
    output_directory = luigi.Parameter()

    def local_target(self, *args):
        return luigi.LocalTarget(
                os.path.join(
                    self.output_directory, *args))
