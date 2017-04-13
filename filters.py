"""
Defines filters for non-linear signals.
"""
from functools import partial
from math import copysign
from operator import add, mul


# Utilities


def multifilter(filters, reductop=all):
    """
    Create a composite filter.
    
    Generate a function that applies multiple filters at once. The combination
    of the results is controlled by the 'reductop' argument.
    Defaults to 'all', ie: all filters must evalute to True.
    
    Parameters
    ----------
    filters : list of callables
        The filters to apply.
    reductop : callable, default `all`
        Reduction operator used to combine the results of individual filters.
        Common choices include `all` (for AND) and `any` (for OR).
    
    Returns
    -------
    composite: callable
        Composite filter.
    """
    filters = tuple(filters)  # make immutable -  grab a copy.
    def composite(k):
        return reductop([f(k) for f in filters])
    return composite


# Index filters


def emit_along(d, k):
    """
    Index filter for phase matching direction k_d.
    
    Parameters
    ----------
    d : int
        Direction to match (ie: `3` for k_3)
    k : list of int
        Interaction indices (ie: `[-1,2,3]` for -k_1+k_2+k_3.
    
    Returns
    -------
    accept : bool
        k is phase-matched.
    """
    return ((sum([copysign(1, i) for i in k if abs(i)==d]) == 1)
            and (sum([copysign(1, i) for i in k if abs(i)!=d]) == 0))


def phasematch_filter(d):
    """
    Make index filter that phase matches wavevector k_d.
    
    Essentially a wrapper over `emit_along`.
    """
    def emit_along_d(k):
        return emit_along(d, k)
    return emit_along_d


def fundamental_freq(k):
    """
    Index filter for fundamental spectral band.
    
    Total interaction order is 1 such that the spectral frequency is close
    to the input.
    
    Parameters
    ----------
    k : list of int
        Interaction indices (ie: `[-1,2,3]` for -k_1+k_2+k_3.
    
    Returns
    -------
    accept : bool
    """
    return abs(sum([copysign(1,i) for i in k])) == 1


_pp_filter = multifilter([fundamental_freq, phasematch_filter(3)])


def pump_probe(k):
    """
    Filter for 3-pulses pump-probe geometry.
    
    It requires phase matching along k_3 and the total interaction order
    to stay 1.
    """
    return _pp_filter(k)


# term filters

def filter_terms(f, expr):
    return expr.func(*filter(f, expr.args))


def contains(s, sig):
    """Signal contains string s"""
    return s in repr(sig)


def contains_filter(s):
    return partial(contains, s)