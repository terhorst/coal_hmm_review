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
        self.configs = [k for k in self.sfs if k is not None]
        self.config_array = momi.config_array(populations, self.configs)

    # def events(self, t_m, t_div, log_p, N_0, N_1, N_A):
    def events(self, t_div, N_0, N_1, N_A):
        events = [
                ('-en', 0., self.populations[0], N_0), 
                ('-en', 0., self.populations[1], N_1), 
                ('-ej', t_div, self.populations[1], self.populations[0]),
                ('-en', t_div, self.populations[0], N_A), 
                ]
        return events
        dm = self.scale(GlobalConfig().migration_density)
        t = t_m + dm
        while t <= t_div:
            for i in (0, 1):
                events.append(('-ep', t, self.populations[i], self.populations[1 - i], 10 ** log_p))
            t += dm
        return events

    def negloglik(self, x):
        # t_m, t_div_inc, N_0, N_1, N_A, log_p = x
        # t_div = t_m + t_div_inc
        # events = self.events(t_m, t_div, log_p, N_0, N_1, N_A)
        events = self.events(*x)
        demo = momi.make_demography(events,
                                    sampled_pops=self.populations,
                                    sampled_n=self.n,
                                    default_N=1,
                                    time_scale="standard")
        try:
            eSFS = momi.expected_sfs(demo, self.config_array, mut_rate=1.0)
            eTot = momi.expected_total_branch_len(demo)
        except AssertionError:
            print("got inf from momi")
            return np.inf
        sdict = {k: v for k, v in zip(self.configs, eSFS)}
        sv = sum(self.sfs.values())
        ll = sum(k * anp.log(sdict[sp]) for sp, k in self.sfs.items())
        ll -= sv * THETA * eTot
        v = []
        for sp in self.sfs:
            k = self.sfs[sp]
            v.append((k, sdict[sp] / eTot, 1. * k / sv, sp))
        # pprint(sorted(v))
        assert ll < 0
        print(x, ll)
        return -ll

    def scale(self, x):
        return x / 29. / (2. * N0)

    def run(self):
        # bounds = [(0., self.scale(3e5)), 
        #           (self.scale(1e2), self.scale(3e5))]
        # bounds += [(1e2 / (2. * N0), 1e7 / (2. * N0))] * 3
        # bounds += [(-6, -1)]
        # x0 = [self.scale(1e4), self.scale(1e5), 10., 10., 10., -3]
        # def f_no_p(x, log_p):
        #     xx = list(x) + [log_p]
        #     return self.negloglik(xx)
        # res = scipy.optimize.minimize(
        #         f_no_p,
        #         jac=autograd.grad(f_no_p, 0), 
        #         bounds=bounds[:-1],
        #         x0=x0[:-1], 
        #         args=(x0[-1],),
        #         method="TNC")
        # print(res)
        # x0[:-1] = res.x
        # def f_p(log_p, x):
        #     print(log_p, x)
        #     xx = list(x) + [log_p]
        #     return self.negloglik(xx)
        # res = scipy.optimize.minimize_scalar(
        #         f_p,
        #         bounds=bounds[-1],
        #         method="bounded",
        #         args=(x0[:-1],))
        # print(res)
        # x0[-1] = res.x
        bounds = [(self.scale(1e3), self.scale(1e6))] + [(.01, 100)] * 3
        x0 = [self.scale(1e5), 1, 1, 1]
        res = scipy.optimize.minimize(self.negloglik, 
                jac=autograd.grad(self.negloglik, 0), 
                bounds=bounds, x0=x0, method="TNC")
        print(res)
        return TwoPopulationDemography(
                params=res.x,
                n=self.n, 
                populations=self.populations, 
                events=self.events(*res.x))
