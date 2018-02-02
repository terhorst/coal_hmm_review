import luigi
import pickle
import sh
import os
import os.path

smc = sh.Command("smc++")

import data.original
import config
import tasks

class Pop2SMC(luigi.Task):
    chromosome = luigi.IntParameter()
    pop1 = luigi.Parameter()
    pop2 = luigi.Parameter(default=None)

    def requires(self):
        return {'pops': tasks.PopulationMap(),
                'vcf': data.original.OriginalFullVCF(chromosome=self.chromosome)}

    def output(self):
        dir = self.pop1
        if self.pop2 is not None:
            dir += "-" + self.pop2
        return config.GlobalConfig().local_target("smc", dir,
                                                  'data', 
                                                  "chr%d.smc.gz" % self.chromosome)

    def run(self):
        self.output().makedirs()
        pops = pickle.load(open(self.input()['pops'].path, 'rb'))
        samples = [(self.pop1, pops[self.pop1])]
        if self.pop2 is not None:
            samples.append((self.pop2, pops[self.pop2]))
        smc.vcf2smc(self.input()['vcf'].path, self.output().path,
                    self.chromosome,
                    *["%s:%s" % (p1, ",".join(p2)) for p1, p2 in samples])


class RunSMCMarginal(luigi.Task):
    pop = luigi.Parameter()
    def requires(self):
        return [Pop2SMC(pop1=self.pop, chromosome=c) for c in range(1, 23)]
    def output(self):
        return config.GlobalConfig().local_target("smc", self.pop, 'estimate', "model.final.json")
    def run(self):
        self.output().makedirs()
        smc.estimate("-o", os.path.split(self.output().path)[0], "-v", 1.25e-8, 
                     "-c", 10000, 
                     *[f.path for f in self.input()],
		     cores=4)


class RunSMCJoint(luigi.Task):
    pop1 = luigi.Parameter()
    pop2 = luigi.Parameter()

    def requires(self):
        keys = [{'pop1': self.pop1}, {'pop1': self.pop2},
                {'pop1': self.pop1, 'pop2': self.pop2},
                {'pop2': self.pop1, 'pop1': self.pop2}]
        return {'pop1': RunSMCMarginal(self.pop1),
                'pop2': RunSMCMarginal(self.pop2),
                'data': [Pop2SMC(chromosome=c, **k)
                         for k in keys
                         for c in range(1, 23)]}

    def output(self):
        return config.GlobalConfig().local_target('smc', 
                                                  '%s-%s' % (self.pop1, self.pop2),
                                                  'split',
                                                  'model.final.json')

    def run(self):
        self.output().makedirs()
        smc.split('-o', os.path.split(self.output().path)[0],
                  self.input()['pop1'].path,
                  self.input()['pop2'].path,
                  *[f.path for f in self.input()['data']])
