"""
Defines non-linear signals and utility functions.
"""
# written for python 2

from __future__ import division, print_function

from math import copysign
from sympy.core.symbol import Symbol
from sympy.core.cache import cacheit
import sympy as sy
import operator
from functools import reduce
from itertools import product, combinations_with_replacement, chain, repeat
try:
    from itertools import ifilter as filter
except ImportError:
    pass
import warnings


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


def phase(i):
    "Phase component of amplitude(i)"
    assert i != 0
    s = sy.symbols("\phi_"+str(abs(i)), real=True)
    if i < 0:
        s = -s
    return s


def norm(i):
    "Norm of amplitude(i)"
    assert i != 0
    s = sy.symbols("a_"+str(abs(i)), real=True)
    return s


def signal(k):
    """
    Generate a complete signal, including phase and power dependence.
    """
    amps = [amplitude(i) for i in k]
    sig = Signal(k)
    return reduce(operator.mul, amps+[sig], 1)


def e_field(i, amp=True):
    """
    Generate the electric field from wavevector index i.
    
    Complex conjugate is returned if i<0.
    """
    ia = abs(i)
    s = sy.symbols("E_"+str(ia))
    if amp:
        s *= amplitude(ia)
    if i < 0:
        s = sy.functions.conjugate(s)
    return s


def strict_ordering(k):
    return all([abs(i)<=abs(j) for i, j in zip(k[:-1], k[1:])])


def signals_for_order(n, n_pulses, strict=True, filters=None):
    """
    Generate signals for orders n for n_pulses. Uses strict ordering by default.

    Returns
    -------
    expr : sympy expression
    """
    # TODO: better doc
    # TODO: support filters is an iterable or a function.
    warnings.warn("Signals for order now doesn't do the sum")
    #f = lambda k, n: list(product(*repeat(k, n)))
    k = range(1, n_pulses + 1)
    k = list(chain(*zip(k, [-i for i in k])))  # [1, -1, 2, -2, ...]
    # take all possible combinations of k with n repeats
    k = list(product(*repeat(k, n))) # maybe product(k, n)
    if strict:
        k = filter(strict_ordering, k)
    if filters:
        k = filter(filters, k)
    return [signal(i) for i in k]


#  expression parsing and manipulation


def list_indices(expr):
    "Gathers the list of pulse indices present in expr. Returns a sorted list"
    return sorted(list(
        set(chain(*[map(abs, a.k)
                    for a in sy.preorder_traversal(expr)
                    if isinstance(a, Signal)]
                  ))))


def max_index(expr):
    """Get largest pulse index"""
    return max(list_indices(expr))


def expand_amps(expr, phase_only=True):
    """Expands the amplitudes to polar coordinates."""
    indices = list_indices(expr)
    expr = expr.subs(
        [(amplitude(i), norm(i)*sy.exp(-sy.I*phase(i)))
         for i in indices])
    if phase_only:
        expr = expr.subs(
            [(norm(i), 1) for i in indices])
    return expr


def weights(alpha, phases):
    """
    Compute phase cycling weights to isolate a given coherence transfer pathway.
    
    Parameters
    ----------
    alpha: (m,) list of ints
        Coherence tranfer pathway. ex: (-1, 1, 1) for rephasing.
    phases: (n,m) list of phases
        Table of phases for n-step phase cycling (rows) with m pulses (columns).
        
    Returns
    -------
    weights: (n,) list of weights (complex)
        List of weights for n-step phase cycling.
    """
    return [sy.exp(-sy.I*sum([a*p
                              for a, p in zip(alpha, pulses)]))
            for pulses in phases]


def coh_transfer(k, n_pulses=None):
    """
    Compute coherence transfer pathway \alpha from interaction indices.
    
    This translates the [-1, 2, 3] convention used here to [-1, 1, 1] used in
    phase cycling theory, where the position indicates the pulse and the value
    indicates the total coherence transfer.
    
    Note that this relation cannot be inverted:
    Both the pump probe k=[-1, 1, 3] and the linear k=[3] yield alpha = [0,0,1]
    """
    if n_pulses is None:
        n_pulses = max([abs(v) for v in k])
    alpha = [0 for i in range(n_pulses)]  # python is weird
    for v in k:
        i = abs(v)-1
        t = copysign(1, v) # keep sign only
        alpha[i] += t
    return alpha


def phase_cycle(expr, phases, weights, norms=None, expand_norms=False):
    """
    Apply phase cycling to expr.

    This applies phase cycling to `expr`. The phases are cycled following the
    table `phases` and optionally the chopping amplitudes `norms`. The resulting
    expressions are summed with the given weights. If weights are specified
    (default), they are determined automatically according to!!!

    Parameters
    ----------
    expr: detected signal
        Expression to phase-cycle.
    phases: (n,m) list of list of phases (reals).
        Table of phases for n-step phase cycling (rows) with m pulses (columns)
    weights: (n,) list of weights (complex)
        List of weights for n-step phase cycling.
    norms: (n,m) list of list of chopping factors (reals)
        Table of amplitude/chopping factors for n-step phase cycling (rows) with
        m pulses (columns). Optional. Defaults to no chopping.

    """
    if norms is None:
        if expand_norms:
            norms = [[norm(i+1) for i in range(len(phases[0]))]
                     for j in range(len(phases))]
        else:
            norms = repeat(repeat(1))
    cycled = sum([w * expr.subs(
                        [(amplitude(i + 1), a * sy.exp(-sy.I * p))
                         for i, (a, p) in enumerate(zip(rads, phis))])
                  for w, rads, phis in zip(weights, norms, phases)])
    return cycled#.expand() # this will simplify things