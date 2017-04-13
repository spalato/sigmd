"""
Defines non-linear signals and utility functions.
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


def e_field(i):
    """
    Generate the electric field from wavevector index i.
    
    Complex conjugate is returned if i<0.
    """
    i = abs(i)
    s = amplitude(i)*sy.symbols("E_"+str(i))
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
    #f = lambda k, n: list(product(*repeat(k, n)))
    k = range(1, n_pulses + 1)
    k = list(chain(*zip(k, [-i for i in k])))  # [1, -1, 2, -2, ...]
    # take all possible combinations of k with n repeats
    k = list(product(*repeat(k, n))) # maybe product(k, n)
    if strict:
        k = filter(strict_ordering, k)
    if filters:
        k = ifilter(filters, k)
    return sum([signal(i) for i in k])


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


def phase_cycle(expr, phases, weights=None, norms=None):
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

    norms: (n,m) list of list of chopping factors (reals)
        Table of amplitude/chopping factors for n-step phase cycling (rows) with
        m pulses (columns). Optional. Defaults to no chopping.

    weights: (n,) list of weights (complex)
        List of weights for n-step phase cycling. By default, they are
        determined automatically.

    """
    # TODO: handle norms, handle weights.
    if norms is None:
        norms = repeat(repeat(1))
    cycled = sum([w * expr.subs(
                        [(amplitude(i + 1), a * sy.exp(-sy.I * p))
                         for i, (a, p) in enumerate(zip(rads, phis))])
                  for w, rads, phis in zip(weights, norms, phases)])
    return cycled