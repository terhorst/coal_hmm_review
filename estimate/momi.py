import numpy as np
import autograd
import autograd.numpy as anp
import scipy.optimize
from pprint import pprint
import momi
import pdb
import attr

from config import *

@attr.s
class TwoPopulationDemography:
    params = attr.ib()
    n = attr.ib()
    populations = attr.ib()
    events = attr.ib()


class PairwiseMomiEstimator:
    """
    Estimate two-population split model with residual migration.
    In this model the free parameters are the divergence time,
    residual migration rate (assumed symmetric), and the stopping
    time of migration.
    """
    def __init__(self, populations, sfs, n):
        self.populations = populations
        self.sfs = {k: v for k, v in sfs.items() if np.array_equal(list(map(sum, k)), n)}
        self.n = n

    # def events(self, t_m, t_div, log_p, N_0, N_1, N_A):
    def events(self, t_div, N_0, N_1, N_A):
        events = [
                ('-en', 0., self.populations[0], N_0), 
                ('-en', 0., self.populations[1], N_1), 
                ('-ej', t_div, self.populations[1], self.populations[0]),
                ('-en', t_div, self.populations[0], N_A), 
                ]
        # dm = self.scale(GlobalConfig().migration_density)
        # t = t_m + dm
        # while t <= t_div:
        #     for i in (0, 1):
        #         events.append(('-ep', t, self.populations[i], self.populations[1 - i], 10 ** log_p))
        #     t += dm
        # print(events)
        return momi.make_demography(
                events,
                sampled_pops=self.populations,
                sampled_n=self.n,
                default_N=1.,
                time_scale="standard")

    def scale(self, x):
        return x / 29. / (2. * N0)

    def run(self):
        # bounds = [(self.scale(1e3), self.scale(1e6))] * 2 + [(.01, 100)] * 3 + [(-6, -1)]
        # x0 = [self.scale(1e4), self.scale(1e5), 1, 1, 1, -4]
        bounds = [(self.scale(1e3), self.scale(1e6))] * 1 + [(.01, 100)] * 3
        x0 = [self.scale(5e4), 1, 1, 1]
        seg_sites = momi.site_freq_spectrum(self.populations, [self.sfs])
        surface = momi.SfsLikelihoodSurface(seg_sites, demo_func=self.events)
        res = surface.find_mle(x0, bounds=bounds)
        return TwoPopulationDemography(
                params=res.x,
                n=self.n, 
                populations=self.populations, 
                events=self.events(*res.x))
