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
        self.sfs = sfs
        self.n = n
        self.configs = [k for k in sfs if k is not None]
        self.config_array = momi.config_array(populations, self.configs)

    def events(self, t_m, t_div, log_p, N_0, N_1, N_A):
        events = [
                ('-en', 0., self.populations[0], N_0), 
                ('-en', 0., self.populations[1], N_1), 
                ('-ej', t_div, self.populations[1], self.populations[0]),
                ('-en', t_div, self.populations[0], N_A), 
                ]
        t = t_m * (1. + 1e-8)
        while t <= t_div:
            for i in (0, 1):
                events.append(('-ep', t, self.populations[i], 
                               self.populations[1 - i], 
                               10 ** log_p))
            t += self.scale(25000.)
        return events

    def negloglik(self, x):
        t_m, t_div_inc, log_p, N_0, N_1, N_A = x
        t_div = t_m + t_div_inc
        events = self.events(t_m, t_div, log_p, N_0, N_1, N_A)
        demo = momi.make_demography(events,
                                    sampled_pops=self.populations,
                                    sampled_n=self.n,
                                    default_N=1.,
                                    time_scale="standard")
        eSFS = momi.expected_sfs(demo, self.config_array)
        eTot = momi.expected_total_branch_len(demo)
        sdict = {k: THETA * v for k, v in zip(self.configs, eSFS)}
        assert sum(sdict.values()) < .1
        sdict[None] = 1. - THETA * eTot
        print(sdict)
        ll = sum(k * anp.log(sdict[sp]) for sp, k in self.sfs.items())
        print(x, ll)
        return -ll

    def scale(self, x):
        return x / 29. / (2. * N0)

    def run(self):
        bounds = [(0., self.scale(1e5)), 
                  (self.scale(1e2), self.scale(1e5)), 
                  (-8, -1)]
        bounds += [(1e2 / (2. * N0), 1e6 / (2. * N0))] * 3
        x0 = np.mean(bounds, axis=1)
        grad = autograd.grad(self.negloglik, 0)
        res = scipy.optimize.minimize(self.negloglik, jac=grad, bounds=bounds, x0=x0)
        return TwoPopulationDemography(
                params=res.x,
                n=self.n, 
                populations=self.populations, 
                events=self.events(*res.x))
