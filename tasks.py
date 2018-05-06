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
psmc = sh.Command("/scratch/psmc/psmc")
psmcplot = sh.Command("/scratch/psmc/utils/psmc_plot.pl")
psmc2hist = sh.Command("/scratch/psmc/utils/psmc2history.pl")
psmc2msmc = sh.Command(os.path.join(basedir, "scripts", "psmc2msmc.R"))
msmc = sh.Command("/scratch/msmc/msmc_1.0.0_linux64bit")
msmcplot = sh.Command(os.path.join(basedir, "scripts", "msmcplot.R"))
multihet = sh.Command("/scratch/msmc/generate_multihetsep.py").bake(_no_err=True)

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
                self.seed,
                self.N,
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
class VCF2MSMC(luigi.Task):
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
            "estimate",
            f"msmc.final.txt")

    def run(self):
        msmc(*[f.path for f in self.input()],
             i=10,
             outFilePrefix=self.output().path[:-len(".final.txt")])

@luigi.util.requires(EstimateSizeHistoryMSMC)
class PlotMSMC(luigi.Task):
    def output(self):
        return self.local_target(
            "msmc",
            "estimate",
            "plot_msmc.pdf")

    def run(self):
        msmcplot(self.output().path, self.input().path)

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
class SMC2PSMC(luigi.Task):
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
        return self.local_target(
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
        return self.local_target(
            "psmc",
            "data",
            "combined.psmcfa")

    def run(self):
        sh.cat(*[f.path for f in self.input()], _out=self.output().path)

@luigi.util.requires(PSMCCombiner)
class EstimateSizeHistoryPSMC(luigi.Task):
    def output(self):
        return self.local_target(
            "psmc",
            "estimates.txt")

    def run(self):
        psmc(self.input().path, _out=self.output().path)


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


@luigi.util.requires(EstimateSizeHistoryPSMC)
class PlotPSMC(luigi.Task):
    def output(self):
        return self.local_target(
            "psmc",
            "plot_psmc.pdf")

    def run(self):
        psmcplot("-g", 29, "-p", self.output().path[:-len(".pdf")], self.input().path)

class EstimateAll(SimulationTask, luigi.WrapperTask):
    def requires(self):
        return [self.clone(cls) for cls in (PlotPSMC,
                                            PlotMSMC,
                                            PlotSMC)]

class ManyEstimates(luigi.Config):
    N = luigi.IntParameter()
    n_replicates = luigi.IntParameter(default=10)

class EstimateManyReplicates(ManyEstimates, luigi.WrapperTask):
    def requires(self):
        return [self.clone(EstimateAll, demography=demo, seed=i)
                for i in range(1, 1 + self.n_replicates)
                for demo in demography.DEMOGRAPHIES]


class CombinePlotsSMC(ManyEstimates):
    demography = luigi.Parameter()
    def requires(self):
        return [self.clone(EstimateSizeHistoryPSMC, 
                           demography=self.demography,
                           seed=i)
                for i in range(1, 1 + self.n_replicates)]

    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.N,
            "smc_combined.pdf")

    def run(self):
        smc.plot('-g', 29, self.output().path, 
                 *[f.path for f in self.input()])

class CombinePlotsMSMC(ManyEstimates):
    demography = luigi.Parameter()
    def requires(self):
        return [self.clone(EstimateSizeHistoryMSMC, 
                           demography=self.demography,
                           seed=i)
                for i in range(1, 1 + self.n_replicates)]

    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.N,
            "msmc_combined.pdf")

    def run(self):
        msmcplot(self.output().path, 
                 *[f.path for f in self.input()])

class CombinePlotsPSMC(ManyEstimates):
    demography = luigi.Parameter()
    def requires(self):
        return [self.clone(EstimateSizeHistoryPSMC, 
                           demography=self.demography,
                           seed=i)
                for i in range(1, 1 + self.n_replicates)]

    def output(self):
        return GlobalConfig().local_target(
            self.demography,
            self.N,
            "psmc_combined.pdf")

    def run(self):
        hists = []
        for f in self.input():
            hists.append(f.path + '.hist')
            psmc2msmc(psmc2hist(f.path), _out=hists[-1])
