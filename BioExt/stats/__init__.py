
from __future__ import division, print_function

import numpy as np

from scipy.misc import comb, factorial
from scipy.optimize import fmin_cobyla
from scipy.stats import rv_continuous, rv_discrete
from scipy.special import beta, gamma


__all__ = ['bbinom', 'pois']


class bbinom_gen(rv_discrete):

    def _pmf(self, k, n, a, b):
        return comb(n, k) * beta(k + a, n - k + b) / beta(a, b)

    def fit(self, data, *args, **kwds):
        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)
        if floc is None and fscale is None:
            if not isinstance(data, np.ndarray):
                data = np.array(data)
            n, a0, b0 = bbinom_gen._fitstart(data)
            tot = np.sum(data)
            # function to minimize
            ks = np.arange(n + 1)
            def f(x):
                return np.sum( (data - tot * bbinom_gen._pmf(self, ks, n, x[0], x[1])) ** 2 )
            # initial guess with method of moments
            x0 = np.array([a0, b0])
            xn = fmin_cobyla(f, x0, [lambda x: x[0], lambda x: x[1]], rhobeg=1e-1, rhoend=1e-10, disp=0)
            return n, xn[0], xn[1]
        else:
            raise ValueError('floc and fscale unsupported')

    @staticmethod
    def _fitstart(data):
        n = len(data) - 1
        tot = sum(data)
        # method of moments initial estimation
        m1 = np.sum(np.arange(n + 1) * data) / tot
        m2 = np.sum((np.arange(n + 1) ** 2) * data) / tot
        ah = (n * m1 - m2) / (n * ((m2 / m1) - m1 - 1) + m1)
        bh = (n - m1) * (n - (m2 / m1)) / (n * ((m2 / m1) - m1 - 1) + m1)
        return n, ah, bh

bbinom = bbinom_gen(name='bbinom', shapes='n, a, b')


class pois_gen(rv_discrete):

    def _pmf(self, k, lam):
        return (lam ** k) / factorial(k) * np.exp(-lam)

    def fit(self, data, *args, **kwds):
        floc = kwds.get('floc', None)
        fstart = kwds.get('fstart', None)
        if floc is None and fstart is None:
            if not isinstance(data, np.ndarray):
                data = np.array(data)
            lam0, = pois_gen._fitstart(data)
            tot = np.sum(data)
            # function to minimize
            ks = np.arange(len(data))
            def f(x):
                return np.sum( (data - tot * pois_gen._pmf(self, ks, x[0])) ** 2 )
            # initial guess with method of moments
            x0 = np.array([lam0])
            xn = fmin_cobyla(f, x0, [lambda x: x[0]], rhobeg=1e-1, rhoend=1e-10, disp=0)
            return (xn[0],)
        else:
            raise ValueError('floc and fscale unsupported')

    @staticmethod
    def _fitstart(data):
        return (np.sum(np.arange(len(data)) * data) / np.sum(data),)

pois = pois_gen(name='pois', shapes='lam')


# class diri_gen(rv_continuous):
#
#     def _pdf(self, x, a):
#         return gamma(np.sum(a)) / np.prod(gamma(a)) * np.prod(np.power(x, a - 1))
#
# diri = diri_gen(name='diri', shapes='a')
