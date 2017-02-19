import luigi
import sh
import os
import pysam
import collections
import pickle

from original_data import VCFTask
import config

smc = sh.Command("smc++")

class ConvertAll(VCFTask):
    'Catch-all task to batch convert everything at once'

    def run(self):
        # This must be dynamic because we depend on self.populations
        tasks = []
        for chrom in range(1, 23):
            chrom = str(chrom)
            tasks.append(VCF2Momi(contig=chrom))
            for pop in self.populations:
                tasks.append(VCF2SMC(contig=chrom, population=pop))
        yield tasks


class VCF2SMC(VCFTask):
    population = luigi.Parameter()
    contig = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(
                os.path.join(
                    config.GlobalConfig().output_directory, "smc", 
                    self.population, self.contig + ".txt.gz"))
        
    def run(self):
        self.output().makedirs()
        samples = self.populations[self.population]
        smc("vcf2smc", "-m", self.input()['centromeres'].path,
                "--drop-first-last",
                self.input()['vcf'].path, 
                self.output().path,
                self.contig,
                "{}:{}".format(self.population, ",".join(samples)))
         

class VCF2Momi(VCFTask):
    contig = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(config.GlobalConfig().output_directory, "momi", self.contig + ".dat"))

    def run(self):
        self.output().makedirs()
        sfs = collections.Counter()
        pd = self.populations
        pops = list(self.populations)
        with pysam.VariantFile(self.input()['vcf'].path) as vcf:
            for record in vcf.fetch(contig=self.contig):
                d = {}
                for pop in pops:
                    gts = [x for sample in pd[pop]
                           for x in record.samples[sample]['GT']
                           if x is not None]
                    n = len(gts)
                    a = sum(gts)
                    d[pop] = (n - a, a)
                k = tuple([d[pop] for pop in pops])
                sfs[k] += 1
            sfs[None] = vcf.header.contigs[self.contig].length - sum(sfs.values())
        pickle.dump({'sfs': sfs, 'pops': pops}, open(self.output().path, "wb"), -1)
