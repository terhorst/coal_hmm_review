import numpy as np
import matplotlib
matplotlib.use("Agg")
import luigi.task
import tabulate
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from config import *
import tasks
import estimate

class SummarizePairwiseMLE(tasks._AllPairs):
    def run(self):
        populations = sorted(self.population_map)
        L = len(populations)
        df = pd.DataFrame(np.zeros([L, L]), index=populations, columns=populations)
        pairwise = yield self.all_pairs_tasks()
        for pair in pairwise:
            print(pair)
            data = unpickle(pairwise[pair])
            tm, tdiv, *N, log_p = data.params
            df[pair[0]][pair[1]] = df[pair[1]][pair[0]] = tdiv * 29 * 20

        mask = np.zeros_like(df, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True

        # Set up the matplotlib figure
        f, ax = plt.subplots(figsize=(11, 9))

        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(220, 10, as_cmap=True)

        # Draw the heatmap with the mask and correct aspect ratio
        sns.heatmap(df, mask=mask, cmap=cmap, vmax=200,
                square=True, linewidths=.5, cbar_kws={"shrink": .5},
                ax=ax)
        self.output().makedirs()
        plt.savefig(self.output().path)

    def output(self):
        return {ext:
                GlobalConfig().local_target("summaries", "pairwise_divergence." + ext)
                for ext in ['txt', pdf']}
