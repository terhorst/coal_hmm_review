import luigi
import os
import pysam
import collections
import pickle
import re

from config import *
import util

## Files provided by our collaborators
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

## Derived data sets
class PopulationMap(luigi.Task):
    """
    Dict populations => samples. Restricted to have only samples / pops
    that are actually in the VCF.
    """
    def requires(self):
        return {"vcf": IndexedVCF(), 
                "populations": OriginalPopulations()}

    def output(self):
        return luigi.LocalTarget(
                os.path.join(
                    GlobalConfig().output_directory,
                    "vcf_population_map.dat"))

    def run(self):
        pops = {}
        with open(self.input()['populations'].path, "rt") as f:
            fields = next(f)[1:].strip().split("\t")
            for line in f:
                # The provided spreadsheet is malformatted. We can
                # restrict to the first len(fields) entries and the
                # last will be the population ID, with Mbororo Fulani
                # truncated to Mbororo.
                record = dict(zip(fields,
                                  re.split(r"\s+", line.strip())[:len(fields)]
                                  ))
                pops.setdefault(record['Population'].replace(
                    " ", "_"), []).append(record['SampleID'])
        with pysam.VariantFile(self.input()['vcf'].path) as vcf:
            vcf_samples = set(vcf.header.samples)
        pops = {pop: vcf_samples & set(samples) for pop, samples in pops.items()}
        pops = {pop: samples for pop, samples in pops.items() if samples}
        pickle.dump(pops, open(self.output().path, "wb"), -1)


class _VCFConverter(luigi.Task):
    def requires(self):
        return {"vcf": IndexedVCF(), 
                "centromeres": IndexedCentromeres(),
                "populations": Populations()}
        return ret

    @property
    def populations(self):
        return pickle.load(open(self.input()['populations'].path, "rb"))


class VCF2SMC(_VCFConverter):
    population = luigi.Parameter()
    contig = luigi.Parameter()
    distinguished = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(
                os.path.join(
                    GlobalConfig().output_directory, "smc", "data",
                    self.population, 
                    "{}.{}.txt.gz".format(self.distinguished, self.contig)))
        
    def run(self):
        # Composite likelihood over first 3 individuals
        self.output().makedirs()
        samples = self.populations[self.population]
        undistinguished = set(samples) - set([self.distinguished])
        smc("vcf2smc", "-m", self.input()['centromeres'].path,
                "--drop-first-last",
                '-d', self.distinguished, self.distinguished,
                self.input()['vcf'].path, 
                self.output().path,
                self.contig,
                "{}:{}".format(self.population, ",".join(undistinguished)))
         

class VCF2Momi(_VCFConverter):
    contig = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(GlobalConfig().output_directory, 
                "momi", 
                "data", self.contig + ".dat"))

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
