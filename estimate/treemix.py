import pysam
import itertools as it
import pickle
import luigi
import sh
import numpy as np
import gzip

import data.original
import config
import tasks


class TreeMixChrom(luigi.Task):
    chromosome = luigi.IntParameter()

    def requires(self):
        return {'pops': tasks.PopulationMap(),
                'data': data.original.OriginalFullVCF(chromosome=self.chromosome)}

    def output(self):
        return config.GlobalConfig().local_target("treemix", "data",
                                                  "treemix_chr%d.txt.gz" % self.chromosome)

    def run(self):
        self.output().makedirs()
        pops = pickle.load(open(self.input()['pops'].path, 'rb'))
        with gzip.open(self.output().path, "wt") as f, \
                pysam.VariantFile(self.input()['data'].path) as vcf:
            print(" ".join(pops), file=f)
            for row in vcf.fetch():
                if row.ref in "ACTG" and len(row.alts) == 1 and row.alts[0] in "ACTG":
                    out = []
                    for p in pops:
                        gts = [x for sample in pops[p]
                               for x in row.samples[sample]['GT']
                               if x is not None]
                        n = len(gts)
                        a = sum(gts)
                        t = (a, n - a)
                        out.append(t)
                    if min(np.max(out, axis=0)) > 0:
                        print(" ".join(["%d,%d" % t for t in out]), file=f)


class TreeMixData(luigi.Task):
    def requires(self):
        return [TreeMixChrom(c) for c in range(1, 23)]

    def output(self):
        return config.GlobalConfig().local_target("treemix", "data", "combined.txt.gz")

    def run(self):
        self.output().makedirs()
        with gzip.open(self.output().path, "wt") as out:
            for i, c in enumerate(self.input()):
                sh.tail(sh.zcat(c.path), "-n", "+" + str(min(1, i)), out_=out)
