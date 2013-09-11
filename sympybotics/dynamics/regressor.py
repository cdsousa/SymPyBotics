from copy import copy, deepcopy
from sympy import zeros, Matrix
from ..utils import identity
from .simple_fric import frictionterm
from .rne import rne_forward, rne_backward


def regressor(rbtdef, geom, ifunc=None):
    '''Generate regression matrix.'''

    if not ifunc:
        ifunc = identity

    fw_results = rne_forward(rbtdef, geom, ifunc)

    rbtdeftmp = deepcopy(rbtdef)

    dynparms = rbtdef.dynparms()

    Y = zeros((rbtdef.dof, len(dynparms)))

    if rbtdef.frictionmodel == 'simple':
        fric = frictionterm(rbtdef)
        fric_dict = dict(zip(rbtdef.fc, [0] * len(rbtdef.fc)))
        fric_dict.update(dict(zip(rbtdef.fv, [0] * len(rbtdef.fv))))

    for p, parm in enumerate(dynparms):

        for i in range(rbtdef.dof):
            rbtdeftmp.Le[i] = map(
                lambda x: 1 if x == parm else 0, rbtdef.Le[i])
            rbtdeftmp.l[i] = Matrix(rbtdef.l[i]).applyfunc(
                lambda x: 1 if x == parm else 0)
            rbtdeftmp.m[i] = 1 if rbtdef.m[i] == parm else 0

        Y[:, p] = rne_backward(rbtdeftmp, geom, fw_results, ifunc=ifunc)

        if rbtdef.frictionmodel == 'simple':
            select = copy(fric_dict)
            select.update({parm: 1})
            Y[:, p] = ifunc(Y[:, p] + fric.subs(select))

    return Y
