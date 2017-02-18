import luigi
import sh
import pysam
import re
import os

tabix = sh.Command("tabix")

class Centromeres(luigi.ExternalTask):
    def output(self):
        return luigi.LocalTarget("/scratch/terhorst/tishkoff/centromeres.bed.gz")

class OriginalVCF(luigi.ExternalTask):
    def output(self):
        return luigi.LocalTarget("/scratch/terhorst/tishkoff/PASS_SNP.vcf.gz")

class PopulationIndex(luigi.ExternalTask):
    def output(self):
        return luigi.LocalTarget("/scratch/terhorst/tishkoff/samples_info")

class IndexVCF(luigi.Task):
    def requires(self):
        return OriginalVCF()

    def output(self):
        return luigi.LocalTarget(self.input().path + ".tbi")

    def run(self):
        tabix("-p", "vcf", self.input().path)


class VCFTask(luigi.Task):
    @property
    def populations(self):
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
        return pops

    def _requires(self):
        yield IndexVCF()
        yield luigi.Task._requires(self)

    def requires(self):
        return {'vcf': OriginalVCF(), 
                'populations': PopulationIndex(),
                'centromeres': Centromeres()}
