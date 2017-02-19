import numpy as np
import autograd
import autograd.numpy as anp
import scipy.optimize
import momi

N0 = 1e4
THETA = 1.25e-8 * 2 * N0

class PairwiseMomiEstimator:
    """
    Estimate two-population split model with residual migration.
    In this model the free parameters are the divergence time,
    residual migration rate (assumed symmetric), and the stopping
    time of migration.
    """
    def __init__(self, sfs, n):
        self.sfs = sfs
        self.n = n
        self.configs = [k for k in sfs if k is not None]
        self.config_array = momi.config_array(["pop0", "pop1"], self.configs)

    def events(self, t_M, t_div, log_p):
        events = [('-ej', 'pop1', 'pop0', t_div)]
        t = t_M * (1. + 1e-8)
        while t <= t_div:
            for i in (0, 1):
                events.append(('-ep', "pop%d" % i, "pop%d" % (1 - i), 10 ** log_p))
            t += self.scale(1000.)
        return events

    def negloglik(self, x):
        t_m, t_div_inc, log_p = x
        t_div = t_m + t_div_inc
        events = self.events(t_m, t_div, log_p)
        demo = momi.make_demography(events,
                                    sampled_pops=["pop0", "pop1"],
                                    sampled_n=self.n,
                                    default_N=10000.,
                                    time_scale="standard")
        eSFS = momi.expected_sfs(demo, self.config_array)
        eTot = momi.expected_total_branch_len(demo)
        sdict = {k: THETA * v for k, v in zip(self.configs, eSFS)}
        sdict[None] = 1. - THETA * eTot
        ll = sum(k * anp.log(sdict[sp]) for sp, k in self.sfs.items())
        return -ll

    def scale(self, x):
        return x / 29. / (2. * N0)

    def run(self):
        bounds = [(0., self.scale(300000))] * 2 + [(-8, -1)]
        x0 = np.mean(bounds, axis=1)
        grad = autograd.grad(self.negloglik, 0)
        print(x0)
        return scipy.optimize.minimize(self.negloglik, jac=grad, bounds=bounds, x0=x0)
