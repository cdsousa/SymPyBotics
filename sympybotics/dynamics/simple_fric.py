""" Simple Coulomb and viscouse friction model """

from sympy import zeros, sign
from ..utils import identity


def frictionterm(rbtdef, ifunc=None):
    '''Generate friction forces (simple Coulomb and viscouse model).'''
    if not ifunc:
        ifunc = identity
    fric = zeros((rbtdef.dof, 1))
    for i in range(rbtdef.dof):
        fric[i] = ifunc(rbtdef.fv[i] * rbtdef.dq[
                        i] + rbtdef.fc[i] * sign(rbtdef.dq[i]))
    return fric
