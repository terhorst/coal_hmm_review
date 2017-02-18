import luigi
import sh
import os
import pysam
import collections
import pickle

from original_data import VCFTask
import config

smc = sh.Command("smc++")

class VCF2SMC(VCFTask):
    population = luigi.Parameter()
    contig = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(
                os.path.join(
                    config.GlobalConfig().output_directory, "smc", 
                    self.population) + self.contig + ".txt.gz")
        
    def run(self):
        samples = self.populations[self.population]
        smc("vcf2smc", "-m", self.input()['centromeres'].path,
                self.input()['vcf'].path, 
                self.output().path,
                self.contig,
                "{}:{}".format(self.population, ",".join(samples)))
         

class VCF2Momi(VCFTask):
    contig = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(config.GlobalConfig().output_directory, "momi") + self.contig + ".dat")

    def run(self):
        c = collections.Counter()
        gd = {".": None, "0": 0, "1": 1}
        def gt2tup(gt):
            print(gt)
            return tuple([gd[g] for g in gt[::2]])
        pops = self.populations
        SFSEntry = collections.namedtuple("SFSEntry", pops)
        with pysam.VariantFile(self.input()['vcf'].path) as vcf:
            for record in vcf.fetch(contig=self.contig):
                d = {}
                for pop in pops:
                    gts = [x for sample in pops[pop]
                           for x in record.samples[sample]['GT']
                           if x is not None]
                    n = len(gts)
                    a = sum(gts)
                    d[pop] = (n - a, a)
                c[SFSEntry(**d)] += 1
        pickle.dump(c, open(self.output().path, "wb"), -1)
