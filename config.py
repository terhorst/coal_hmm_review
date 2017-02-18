import luigi

class GlobalConfig(luigi.Config):
    output_directory = luigi.Parameter()
