#!python2

"""
Defines non-linear (and linear) signals.
"""

from __future__ import division, print_function
from sympy.core.symbol import Symbol
from sympy.core.cache import cacheit
import sympy as sy
import operator
from itertools import product, combinations_with_replacement, chain, repeat, ifilter


class Signal(Symbol):
    """
    Signal field, in the sense of non-linear response theory.

    This is a specialization of the basic sympy Symbol class.

    Parameters
    ----------
    k : list of ints
        List of wavevector interactions. A third order response would have a
        list of 3 number. The number is the pulse index and the sign designates
        the direction, ie: <0 implies complex conjugate.
        Thus:
            [-1, 2, 3] is the signal due to -k_1+k_2+k_3, the photon echo
            [-1, 1, 3] is the signal due to -k_1+k_1+k_3, a pump-probe signal
            [1] is the linear signal due to the first pulse, on the ket side.

    Note
    ----
    This is designed to be used by the `signal` function.
    """
    # with help from:
    # https://groups.google.com/forum/#!topic/sympy/pU81Trc_Xr8
    def __new_stage2__(cls, name, k, *a, **kw):
        obj = super(Signal, cls).__xnew__(cls, name, *a, **kw)
        obj.k = k
        return obj

    def __new__(cls, k, *a, **kw):
        s = Signal._mk_symbol(k)
        obj = Signal.__xnew_cached_(cls, s, k, *a, **kw)
        return obj

    __xnew__ = staticmethod(__new_stage2__)
    __xnew_cached_ = staticmethod(cacheit(__new_stage2__))
    # _hashable_content should be fine.
    # (The values of k are included in the symbol.)

    @property
    def order(self):
        return len(self.k)

    @staticmethod
    def _mk_symbol(k):
        """
        Make the string to be used as a symbol, eg: '\chi_{k_1+k_2-k_3}'
        """
        lbl = ['+' if i>0 else '-' for i in k]
        lbl = [s+'k_'+str(abs(n)) for s, n in zip(lbl, k)]
        return r'\chi_{'+''.join(lbl)+'}'

    # yeah, that's more confusing than anything.
    #def _eval_conjugate(self):
        #return -1*Signal([-i for i in self.k])
        # TODO: implement custom complex conjugate: \chi_{k_1} -> \chi_{-k_1}?


def amplitude(i):
    """
    Generate an amplitude from wavevector index i.

    This will be: A_i if i>0, A*_i otherwise.
    """
    assert i != 0
    s = sy.symbols("A_"+str(abs(i)))
    if i > 0:
        return s
    else:
        return sy.functions.conjugate(s)


def signal(k):
    """
    Generate a complete signal, including phase and power dependence.
    """
    amps = [amplitude(i) for i in k]
    sig = Signal(k)
    return reduce(operator.mul, amps+[sig], 1)


def signals_for_order(n, n_pulses, strict=True, filter_=None):
    """
    Generate signals for orders n for n_pulses. Uses strict ordering by default.

    Returns
    -------
    expr : sympy expression
    """
    # TODO: better doc
    # !!!
    # TODO: strict ordering as implemented here may cause a problem.
    # [1, -1] for order 2 can't result in -k1+k1 (only k1+k1, k1-k1, -k1-k1
    # is that a problem? Exemple case: 3rd order 3 pulse: [k3-k3+k3] (SPM)
    # This would influence the creation of amplitudes and
    # responses: they would not commute anymore (maybe?)
    if strict:
        f = combinations_with_replacement
    else:
        f = lambda k, n: list(product(*repeat(k, n)))
    k = range(1, n_pulses+1)
    k = list(chain(*zip(k, [-i for i in k])))  # [1, -1, 2, -2, ...]
    k = f(k, n)
    if filter:
        k = ifilter(filter_, k)
    return sum([signal(i) for i in k])

