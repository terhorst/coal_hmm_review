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
        for v in self.output().values():
            v.makedirs()
        populations = sorted(self.population_map)
        L = len(populations)
        df = pd.DataFrame(np.zeros([L, L]), index=populations, columns=populations)
        pairwise = yield self.all_pairs_tasks()
        with open(self.output()['txt'].path, "wt") as txt:
            print("pop1\tpop2\tt_div", file=txt)
            for pair in pairwise:
                print(pair)
                data = unpickle(pairwise[pair])
                tdiv, *N = data.params
                td = tdiv * 29 * 20
                df[pair[0]][pair[1]] = df[pair[1]][pair[0]] = td
                print("{}\t{}\t{}".format(*pair, td * 1e3), file=txt)
        mask = np.zeros_like(df, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True

        # Set up the matplotlib figure
        f, ax = plt.subplots(figsize=(11, 9))

        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(220, 10, as_cmap=True)

        # Draw the heatmap with the mask and correct aspect ratio
        sns.heatmap(df, mask=mask, cmap=cmap,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5},
                    ax=ax)
        plt.savefig(self.output()['pdf'].path)

    def output(self):
        return {ext:
                GlobalConfig().local_target("summaries", "pairwise_divergence." + ext)
                for ext in ['txt', 'pdf']}
