import luigi
import os.path
import sh
import re
import pickle
import shutil

from config import *

class OriginalFile(luigi.ExternalTask):
    def output(self):
        return luigi.LocalTarget(
                os.path.join(GlobalConfig().input_directory,
                    self.filename))

class Tabixed(luigi.Task):
    """
    Class that ensures that the BGZIPped external dependency also has
    tabix index
    """
    target = luigi.TaskParameter()

    @property
    def filename(self):
        return os.path.basename(self.target.output().path)

    def requires(self):
        return self.target

    def run(self):
        shutil.copyfile(self.target.output().path, self.output().path)
        tabix(self.output().path)

    def complete(self):
        if not os.path.exists(self.output().path + ".tbi"):
            return False
        return luigi.Task.complete(self)

    def output(self):
        return GlobalConfig().local_target(self.filename)
