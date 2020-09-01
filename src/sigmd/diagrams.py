"""
Defines Feynman diagrams and utility functions.
"""
# written for python 3

from .signals import Signal
from sympy.core.cache import cacheit

class SingleSidedFeynmanDiagram(Signal):
    """
    Single-sided feynman diagram.

    Describes the evolution of a coefficient in Hilbert space. To be used in
    composing double-sided feynman diagrams.

    Parameters
    ----------
    k : (n,) ints
        List of wavevector interactions. A third order response would have a
        list of 3 number. The number is the pulse index and the sign designates
        the direction, ie: <0 implies complex conjugate.
        Thus:
            [1, -2, 3] is the signal due to k_1-k_2+k_3
            [1] is the linear signal due to the first pulse, on the ket side.

    c : (n+1,) ints
        List of basis states the diagram goes through. Include initial ground
        state.
    """
    def __new_stage2__(cls, name, k, c, *a, **kw):
        obj = super()._xnew_(cls, name, k, *a, **kw)
        obj.c = c
        return obj

    def __new__(cls, k, c, *a, **kw):
        s = SingleSidedFeynmanDiagram._mk_symbol(k, c)
        obj = Signal.__xnew_cached_(cls, s, k, c, *a, **kw)
        return obj

    __xnew__ = staticmethod(__new_stage2__)
    __xnew_cached_ =  staticmethod(cacheit(__new_stage2__))

    @staticmethod
    def _mk_symbol(k):
        """
        Make the string to be used as a symbol.
        """