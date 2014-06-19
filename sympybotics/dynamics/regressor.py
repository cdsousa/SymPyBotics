from copy import copy, deepcopy
from sympy import zeros, Matrix
from ..utils import identity
from .rne import rne_forward, rne_backward


def regressor(rbtdef, geom, ifunc=None):
    '''Generate regression matrix.'''

    if not ifunc:
        ifunc = identity

    fw_results = rne_forward(rbtdef, geom, ifunc)

    rbtdeftmp = deepcopy(rbtdef)

    dynparms = rbtdef.dynparms()

    Y = zeros(rbtdef.dof, len(dynparms))

    for p, parm in enumerate(dynparms):

        for i in range(rbtdef.dof):
            rbtdeftmp.Le[i] = list(map(
                lambda x: 1 if x == parm else 0, rbtdef.Le[i]))
            rbtdeftmp.l[i] = Matrix(rbtdef.l[i]).applyfunc(
                lambda x: 1 if x == parm else 0)
            rbtdeftmp.m[i] = 1 if rbtdef.m[i] == parm else 0
            rbtdeftmp.Ia[i] = 1 if rbtdef.Ia[i] == parm else 0
            rbtdeftmp.fv[i] = 1 if rbtdef.fv[i] == parm else 0
            rbtdeftmp.fc[i] = 1 if rbtdef.fc[i] == parm else 0
            rbtdeftmp.fo[i] = 1 if rbtdef.fo[i] == parm else 0

        Y[:, p] = rne_backward(rbtdeftmp, geom, fw_results, ifunc=ifunc)

    return Y
