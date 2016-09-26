#!python2

from __future__ import division, print_function
import sympy as sy

a = sy.symbols('a')

ap = sy.conjugate(a)

ap*a
