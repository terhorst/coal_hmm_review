import luigi

import sh
smc = sh.Command("smc++")

class GlobalConfig(luigi.Config):
    input_directory = luigi.Parameter()
    output_directory = luigi.Parameter()
