import pickle
import luigi
import msprime
import contextlib
import tempfile
import collections

from config import *
import data.original
import estimate.tasks
import simulate.msprime

class SimulateTwoPopulationSplitModel(luigi.Task):
    seed = luigi.IntParameter()
    populations = luigi.ListParameter()

    def requires(self):
        return data.original.IndexedVCF()
   
    def output(self):
        return GlobalConfig().local_target("simulation", 
                "-".join(sorted(self.populations)),
                "msprime." + str(self.seed) + ".dat")

    def run(self):
        self.output().makedirs()
        data_set = yield data.tasks.VCF2Momi(vcf_path=self.input().path)
        mle = yield estimate.tasks.PairwiseMomiAnalysis(
                populations=self.populations,
                data_set_paths=[data_set.path])
        d = pickle.load(open(mle.path, "rb"))
        sim = simulate.msprime.MsprimeMomiSimulator(d['populations'], d['n'], d['mle_events'])
        sim.run(self.seed, self.output().path)


class MsprimeToVcf(luigi.Task):
    msprime_dump_path = luigi.Parameter()
    sample_names = luigi.ListParameter()

    def complete(self):
        if not os.path.exists(self.msprime_dump_path + ".vcf.gz.tbi"):
            return False
        return luigi.Task.complete(self)

    def output(self):
        return GlobalConfig().local_target(self.msprime_dump_path + ".vcf.gz")

    def run(self):
        new_samples = ['msp_%d %s' % t for t in enumerate(self.sample_names)]
        tree_seq = msprime.load(self.msprime_dump_path)
        vcf_path = self.output().path[:-3]
        with open(vcf_path, "wt") as f:  # omit the .gz
            tree_seq.write_vcf(f, 2)
        with contextlib.ExitStack() as stack:
            sample_renames = stack.enter_context(tempfile.NamedTemporaryFile("wt"))
            out_vcf_gz = stack.enter_context(open(self.output().path, 'wb'))
            open(sample_renames.name, "wt").write("\n".join(new_samples))
            bgzip(bcftools('reheader', '-s', 
                sample_renames.name, vcf_path, _piped=True), 
                    "-c", _out=out_vcf_gz)
        tabix(self.output().path)
