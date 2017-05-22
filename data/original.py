'Files provided by our collaborators'

import luigi

import util
from config import *

class OriginalVCF(util.OriginalFile):
    filename = "PASS_SNP.vcf.gz"

class OriginalCentromeres(util.OriginalFile):
    filename = "centromeres.bed.gz"

class OriginalPopulations(util.OriginalFile):
    filename = "samples_info"

## Transformed files, with indexes
class IndexedVCF(util.Tabixed):
    target = OriginalVCF()

class IndexedCentromeres(util.Tabixed):
    target = OriginalCentromeres()

class OriginalFullVCF(luigi.ExternalTask):
    chromosome = luigi.IntParameter()

    def run(self):
        tabix(self.output().path)

    def complete(self):
        if not os.path.exists(self.output().path + ".tbi"):
            return False
        return luigi.Task.complete(self)

    def output(self):
        return luigi.LocalTarget(os.path.join(
                GlobalConfig().input_directory,
                'original_data', 'all',
                'all.{}.allsites.vcf.gz'.format(self.chromosome)))
