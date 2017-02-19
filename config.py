import luigi

import sh
smc = sh.Command("smc++")

class GlobalConfig(luigi.Config):
    contigs = list(map(str, range(1, 2)))
    input_directory = luigi.Parameter()
    output_directory = luigi.Parameter()
    
