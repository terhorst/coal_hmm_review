import luigi

import sh
smc = sh.Command("smc++")

class GlobalConfig(luigi.Config):
    output_directory = luigi.Parameter()
