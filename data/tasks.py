import luigi
import luigi.util
import os
import pysam
import collections
import pickle
import re

from config import *
import simulate.tasks
import data.original
import util

## Derived data sets
class PopulationMap(luigi.Task):
    """
    Dict populations => samples. Restricted to have only samples / pops
    that are actually in the VCF.
    """
    def requires(self):
        return {"vcf": data.original.IndexedVCF(), 
                "populations": data.original.OriginalPopulations()}

    def output(self):
        return GlobalConfig().local_target("vcf_population_map.dat")

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
        pops = {pop: list(vcf_samples & set(samples)) for pop, samples in pops.items()}
        pops = {pop: samples for pop, samples in pops.items() if samples}
        pickle.dump(pops, open(self.output().path, "wb"), -1)


class _VCFConverter(luigi.Task):
    def requires(self):
        return {"vcf": data.original.IndexedVCF(), 
                "centromeres": data.original.IndexedCentromeres(),
                "populations": PopulationMap()}
        return ret

    @property
    def populations(self):
        return pickle.load(open(self.input()['populations'].path, "rb"))


class VCF2SMC(_VCFConverter):
    population = luigi.Parameter()
    contig = luigi.Parameter()
    distinguished = luigi.Parameter()

    def output(self):
        return GlobalConfig().local_target(
                    "smc", "data",
                    self.population, 
                    "{}.{}.txt.gz".format(self.distinguished, self.contig))
        
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
                "{}:{}".format(self.population, ",".join(samples)))
         
class _VCF2Momi(_VCFConverter):
    def output(self):
        return luigi.LocalTarget(self.input().path + ".momi.dat")

    def run(self):
        self.output().makedirs()
        sfs = collections.Counter()
        pd = self.populations
        pops = list(self.populations)
        with pysam.VariantFile(self.input().path) as vcf:
            for record in vcf.fetch():
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
        n = {pop: 2 * len(pd[pop]) for pop in pops}
        pickle.dump({'sfs': sfs, 'populations': pops, 'n': n}, 
                     open(self.output().path, "wb"), -1)

class OriginalVCFToMomi(_VCF2Momi):
    def requires(self):
        return data.original.IndexedVCF()

@luigi.util.requires(simulate.tasks.MsprimeToVcf)
class SimulatedVCFToMomi(_VCF2Momi):
    pass
