'Files provided by our collaborators'

import util

class OriginalVCF(util.OriginalFile):
    filename = "PASS_SNP.vcf.gz"

class OriginalCentromeres(util.OriginalFile):
    filename = "centromeres.bed.gz"

class OriginalPopulations(util.OriginalFile):
    filename = "populations.index"

## Transformed files, with indexes
class IndexedVCF(util.Tabixed):
    target = OriginalVCF()

class IndexedCentromeres(util.Tabixed):
    target = OriginalCentromeres()

