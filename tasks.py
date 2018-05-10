import collections
import contextlib
import hashlib
import itertools
import luigi
import luigi.util
import msprime
import numpy as np
import os
import pickle
import pysam
import random
import re
import tempfile

basedir = os.path.dirname(os.path.realpath(__file__))

import sh
smc = sh.Command("smc++")
tabix = sh.Command("tabix")
bgzip = sh.Command("bgzip")
bcftools = sh.Command("bcftools")
psmc = sh.Command("/scratch/psmc/psmc").bake("-N", 15)
psmcplot = sh.Command("/scratch/psmc/utils/psmc_plot.pl").bake("-n", 0) # use best fitting mdl
psmc2hist = sh.Command("/scratch/psmc/utils/psmc2history.pl")
psmc2msmc = sh.Command(os.path.join(basedir, "scripts", "psmc2msmc.R"))
msmc = sh.Command("/scratch/msmc/msmc_1.0.0_linux64bit").bake('-t', 2)
multihet = sh.Command("/scratch/msmc/generate_multihetsep.py").bake(_no_err=True)
msmc2csv = sh.Command(os.path.join(basedir, "scripts", "msmc2csv.R"))
psmc2csv = sh.Command(os.path.join(basedir, "scripts", "psmc2csv.R"))
smc2csv = sh.Command(os.path.join(basedir, "scripts", "smc2csv.R"))
combine_plots = sh.Command(os.path.join(basedir, "scripts", "combine_plots.R")).bake(_no_err=True)

from config import *
import demography
import util
import dical

# All luigi tasks have to go in one big file for now due to circular
# dependencies and luigi.util.require/inherit.

### DATA RELATED TASKS
class SimulationTask(luigi.Config):  
    demography = luigi.Parameter()
    N = luigi.IntParameter()
    seed = luigi.IntParameter()

    @property
    def demo(self):
        return demography.factory(self.demography, self.N)

    def local_target(self, *args):
        return GlobalConfig().local_target(
                self.demography,
                self.N,
                self.seed,
                *args)


class MsPrimeSimulator(SimulationTask):
    contig_id = luigi.Parameter()

    def output(self):
        return self.local_target(f"msp_chr{self.contig_id}.bcf.gz")

    def run(self):
        pc = self.demo.population_configs()
        ev = self.demo.events
        if False:
            msprime.DemographyDebugger(Ne=N0,
                    population_configurations=pc,
                    demographic_events=ev).print_history()
        sim = msprime.simulate(
                mutation_rate=MUTATION_RATE,
                recombination_rate=RECOMBINATION_RATE,
                demographic_events=ev,
                population_configurations=pc,
                random_seed=self.seed + int(self.contig_id) + 1,
                Ne=N0,
                length=GlobalConfig().chromosome_length)
        base = self.output().path[:-len(".bcf.gz")] + ".vcf"
        sim.write_vcf(open(base, "wt"), ploidy=2,
                      contig_id=self.contig_id)
        bcftools.view('-o', self.output().path, '-O', 'b', base)
        bcftools.index('-f', self.output().path)


@luigi.util.requires(MsPrimeSimulator)
class SplitVCF(SimulationTask):
    sample_id = luigi.IntParameter()
    def output(self):
        return self.local_target(
            "msmc",
            "data",
            f"{self.contig_id}_{self.sample_id}.vcf.gz")

    def run(self):
        bcftools.view("-O", "z", "-o", self.output().path, "-s", 
                      f'msp_{self.sample_id}',
                      self.input().path)

@luigi.util.requires(MsPrimeSimulator)
class BCF2VCF(SimulationTask):
    def output(self):
        return self.local_target(f"msp_chr{self.contig_id}.vcf")

    def run(self):
        bcftools.view("-O", "v", "-o", self.output().path, self.input().path)


@luigi.util.inherits(MsPrimeSimulator)
class VCF2MSMC(SimulationTask):
    def requires(self):
        return [self.clone(SplitVCF, sample_id=i) for i in range(2)]

    def output(self):
        return self.local_target(
            "msmc",
            "data",
            f"{self.contig_id}.msmc.txt")

    def run(self):
        multihet(*[f.path for f in self.input()], _out=self.output().path)


class EstimateSizeHistoryMSMC(SimulationTask):
    def requires(self):
        return [self.clone(VCF2MSMC, contig_id=str(k))
                for k in range(GlobalConfig().n_contigs)]

    def output(self):
        return self.local_target(
            "msmc",
            "estimates",
            f"msmc.final.txt")

    def run(self):
        msmc(*[f.path for f in self.input()],
             i=10,
             outFilePrefix=self.output().path[:-len(".final.txt")])

@luigi.util.requires(EstimateSizeHistoryMSMC)
class PlotMSMC(SimulationTask):
    def output(self):
        return self.local_target(
            "msmc",
            "estimates",
            "plot.csv")

    def run(self):
        msmc2csv(MUTATION_RATE, 
                 GENERATION_TIME,
                 self.input().path,
                 self.output().path)

@luigi.util.requires(MsPrimeSimulator)
class VCF2SMC(SimulationTask):
    distinguished = luigi.Parameter()

    def output(self):
        return self.local_target(
            "smc",
            "data",
            f"{self.distinguished}.{self.contig_id}.txt.gz")

    def run(self):
        # Composite likelihood over first 3 individuals
        demo = demography.factory(self.demography, self.N)
        samples = demo.samples()[0] # assume 1 pop for now
        undistinguished = set(samples) - set([self.distinguished])
        self.output().makedirs()
        smc("vcf2smc",
                # "-m", self.input()['centromeres'].path,
                '-d', self.distinguished, self.distinguished,
                self.input().path,
                self.output().path,
                self.contig_id,
                "{}:{}".format('pop1', ",".join(samples)))  # FIXME: hardcoded pop for now

@luigi.util.requires(VCF2SMC)
class SMC2PSMC(SimulationTask):
    def output(self):
        return self.local_target(
            "psmc",
            "data",
            f"{self.distinguished}.{self.contig_id}.psmcfa")

    def run(self):
        data = np.loadtxt(self.input().path, dtype=int)
        L = data[:, 0].sum()
        L += 100 - (L % 100)
        fa = np.full(L, -1)
        last = 0
        for span, a, b, nb in data:
            fa[last:last + span] = a
            last += span
        fa.shape = (L // 100, -1)
        code = fa.max(axis=1).astype('|S1')
        code[code == b'0'] = b'T'
        code[code == b'1'] = b'K'
        code[code == b'2'] = b'T'  # recode monomorphic sites
        code[fa.min(axis=1) == -1] = b'N'
        with open(self.output().path, "wt") as f:
            def fp(x): print(x, file=f)
            fp(">" + self.contig_id)
            Lp = len(code) // 79
            if Lp > 0:
                out = np.full([Lp, 80], b"\n", dtype='string_')
                out[:, :-1] = code[:(79 * Lp)].reshape(Lp, 79)
                fp(out.tostring().decode('ascii')[:-1])  # omit trailing newline
            fp(code[(79 * Lp):].tostring().decode('ascii'))

class EstimateSizeHistorySMC(SimulationTask):
    def requires(self):
        return [self.clone(VCF2SMC, contig_id=str(k),
                           distinguished="msp_0")
                for k in range(GlobalConfig().n_contigs)]

    def output(self):
        return self.local_target(
            "smc",
            "estimates",
            "model.final.json")

    def run(self):
        # Create data sets from composite likelihood
        samples = self.demo.samples()
        out_path = os.path.dirname(self.output().path)
        initial_path = os.path.join(out_path, "initial")
        smc.estimate('-v', 
            '-o', initial_path,
            '--em-iterations', 2,
            '--knots', 24,
            "--cores", 2,
            1.25e-8,
            *[f.path for f in self.input()])
        smc.estimate('-v',
            '-o', out_path,
            '--initial-model', os.path.join(initial_path, 'model.final.json'),
            "--cores", 2,
            "--timepoints", "200,100000",
            1.25e-8,
            *[f.path for f in self.input()])

@luigi.util.requires(EstimateSizeHistorySMC)
class PlotSMC(SimulationTask):
    def output(self):
        return self.local_target(
            "smc",
            "estimates",
            "plot.csv")

    def run(self):
        pdf = self.output().path[:-3] + "pdf"
        smc.plot('-g', GENERATION_TIME, "-c", 
                 pdf,
                 self.input().path)
        smc2csv(self.output().path, 
                self.output().path)

class PSMCCombiner(SimulationTask):
    def requires(self):
        return [self.clone(SMC2PSMC, contig_id=str(k), 
                           distinguished="msp_0")
                for k in range(GlobalConfig().n_contigs)]

    def output(self):
        return self.local_target(
            "psmc",
            "data",
            "combined.psmcfa")

    def run(self):
        sh.cat(*[f.path for f in self.input()], _out=self.output().path)

@luigi.util.requires(PSMCCombiner)
class EstimateSizeHistoryPSMC(SimulationTask):
    def output(self):
        return self.local_target(
            "psmc",
            "estimates.txt")

    def run(self):
        psmc("-o", self.output().path, self.input().path)


class DicalRef(SimulationTask):
    def output(self):
        return self.local_target(
            "dical",
            "ref.fa")

    def run(self):
        L = GlobalConfig().chromosome_length
        open(self.output().path, "wt").write(
            "> contig1\n" +
            ("A" * 80 + "\n") * (L // 80) +
            ("T" * (L % 80))
        )

class EstimateSizeHistoryDical(SimulationTask):
    def requires(self):
        return {
                'vcf': [self.clone(BCF2VCF, contig_id=str(k))
                        for k in range(GlobalConfig().n_contigs)],
                'ref': self.clone(DicalRef)
                }

    def output(self):
        return self.local_target(
            "dical",
            "analysis.out")

    def run(self):
        basename = self.output().path[:-len(".out")]
        da = dical.PieceWiseConstantAnalysis(
                uniqueBasename=basename,
                numCores=2,
                vcfFiles=",".join([f.path for f in self.input()['vcf']]),
                refFiles=",".join([self.input()['ref'].path 
                                  for _ in range(len(self.input()['vcf']))]),
                sampleSize=self.N,
                randomSeed=self.seed
        )
        cmd = da.run()
        os.system(cmd)
        point, lik = da.returnMLE()
        print(point, file=open(basename + ".mle", "wt"))


@luigi.util.requires(EstimateSizeHistoryDical)
class PlotDical(SimulationTask):
    def output(self):
        return self.local_target(
            "dical",
            "plot.csv")

    def run(self):
        # mattias: please fill in
        pass


@luigi.util.requires(EstimateSizeHistoryPSMC)
class PlotPSMC(SimulationTask):
    def output(self):
        return self.local_target(
            "psmc",
            "plot.csv")

    def run(self):
        base = self.output().path[:-len(".csv")]
        psmcplot("-g", GENERATION_TIME, 
                 "-R", 
                 "-u", MUTATION_RATE,
                 "-p", 
                 base,
                 self.input().path)
        psmc2csv(base + ".0.txt", self.output().path)


class PlotCombined(luigi.Task):
    demography = luigi.Parameter()
    N = luigi.IntParameter()
    n_replicates = luigi.IntParameter()

    def requires(self):
        return [self.clone(cls, demography=self.demography, seed=1 + i)
                for i in range(self.n_replicates)
                for cls in (PlotPSMC, PlotMSMC, PlotSMC, PlotDical)[:-1]]

    def output(self):
        return GlobalConfig().local_target(
                self.demography,
                self.N,
                f"{self.demography}.pdf")

    def run(self):
        demo = demography.factory(self.demography, self.N)
        truth_csv = self.output().path[:-4] + "_truth.csv"
        open(truth_csv, "wt").write(demo.to_csv(GENERATION_TIME))
        combine_plots(self.output().path,
                      truth_csv,
                      *[f.path for f in self.input()])


class ManyEstimates(luigi.Config):
    N = luigi.IntParameter()
    n_replicates = luigi.IntParameter(default=10)


class EstimateManyReplicates(ManyEstimates, luigi.WrapperTask):
    def requires(self):
        return [self.clone(EstimateAll, demography=demo, seed=i)
                for i in range(1, 1 + self.n_replicates)
                for demo in demography.DEMOGRAPHIES]
