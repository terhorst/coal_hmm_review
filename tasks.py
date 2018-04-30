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

import sh
smc = sh.Command("smc++")
tabix = sh.Command("tabix")
bgzip = sh.Command("bgzip")
bcftools = sh.Command("bcftools")
psmc = sh.Command("/scratch/psmc/psmc")
psmcplot = sh.Command("/scratch/psmc/utils/psmc_plot.pl")
msmc = sh.Command("/scratch/msmc/msmc_1.0.0_linux64bit")
msmcplot = sh.Command("/scratch/msmc/plot.R")
multihet = sh.Command("/scratch/msmc/generate_multihetsep.py").bake(_no_err=True)

from config import *
import demography
import util

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

class MsPrimeSimulator(SimulationTask):
    contig_id = luigi.Parameter()

    def output(self):
        return GlobalConfig().local_target(
                self.demography,
                self.seed,
                self.N,
                f"msp_chr{self.contig_id}.bcf.gz")

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
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "msmc",
            "data",
            f"{self.contig_id}_{self.sample_id}.vcf.gz")

    def run(self):
        bcftools.view("-O", "z", "-o", self.output().path, "-s", 
                      f'msp_{self.sample_id}',
                      self.input().path)

@luigi.util.inherits(MsPrimeSimulator)
class VCF2MSMC(luigi.Task):
    def requires(self):
        return [self.clone(SplitVCF, sample_id=i) for i in range(2)]

    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
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
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "msmc",
            "estimate",
            f"msmc.final.txt")

    def run(self):
        msmc(*[f.path for f in self.input()],
             i=10,
             outFilePrefix=self.output().path[:-len(".final.txt")])

@luigi.util.requires(EstimateSizeHistoryMSMC)
class PlotMSMC(luigi.Task):
    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "msmc",
            "estimate",
            "plot_msmc.pdf")

    def run(self):
        msmcplot(self.input().path, self.output().path)

@luigi.util.requires(MsPrimeSimulator)
class VCF2SMC(SimulationTask):
    distinguished = luigi.Parameter()

    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
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
class SMC2PSMC(luigi.Task):
    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
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
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "smc",
            "estimates",
            "model.final.json")

    def run(self):
        # Create data sets from composite likelihood
        samples = self.demo.samples()
        smc.estimate('-v', '-o', 
            os.path.dirname(self.output().path),
            "--cores", 2,
            "--knots", 12,
            "--timepoints", "33,100000",
            "--spline", "cubic",
            "-rp", 5,
            1.25e-8,
            *[f.path for f in self.input()])

@luigi.util.requires(EstimateSizeHistorySMC)
class PlotSMC(SimulationTask):
    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "smc",
            "estimates",
            "plot_smc.pdf")

    def run(self):
        smc.plot('-g', 29, self.output().path, self.input().path)

class PSMCCombiner(SimulationTask):
    def requires(self):
        return [self.clone(SMC2PSMC, contig_id=str(k), 
                           distinguished="msp_0")
                for k in range(GlobalConfig().n_contigs)]

    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "psmc",
            "data",
            "combined.psmcfa")

    def run(self):
        sh.cat(*[f.path for f in self.input()], _out=self.output().path)

@luigi.util.requires(PSMCCombiner)
class EstimateSizeHistoryPSMC(luigi.Task):
    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "psmc",
            "estimates.txt")

    def run(self):
        psmc(self.input().path, _out=self.output().path)

@luigi.util.requires(EstimateSizeHistoryPSMC)
class PlotPSMC(luigi.Task):
    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.seed,
            self.N,
            "psmc",
            "plot_psmc.pdf")

    def run(self):
        psmcplot("-p", self.output().path[:-len(".pdf")], self.input().path)

class EstimateAll(SimulationTask, luigi.WrapperTask):
    def requires(self):
        return [self.clone(cls) for cls in (PlotPSMC,
                                            PlotMSMC,
                                            PlotSMC)]

class EstimateManyReplicates(luigi.WrapperTask):
    N = luigi.IntParameter()
    n_replicates = luigi.IntParameter(default=10)

    def requires(self):
        return [self.clone(EstimateAll, demography=demo, seed=i)
                for i in range(1, 1 + self.n_replicates)
                for demo in demography.DEMOGRAPHIES]
