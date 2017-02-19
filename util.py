import luigi
import os.path
import sh
import re
import pickle

from config import *

tabix = sh.Command("tabix")

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
        return self.target.output().path

    def requires(self):
        return self.target

    def run(self):
        tabix(self.filename)

    def complete(self):
        if not os.path.exists(self.filename + ".tbi"):
            return False
        return luigi.Task.complete(self)

    def output(self):
        return luigi.LocalTarget(
                os.path.join(GlobalConfig().input_directory,
                    self.filename))
