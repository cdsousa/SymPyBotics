from copy import deepcopy
from sympy import zeros
from ..utils import identity
from ..geometry import Geometry
from .rne_park import rne_park_forward, rne_park_backward
from .rne_khalil import rne_khalil_forward, rne_khalil_backward


def rne_forward(rbtdef, geom, ifunc=None):
    if rbtdef._dh_convention == 'standard':
        rne_forward = rne_park_forward
    elif rbtdef._dh_convention == 'modified':
        rne_forward = rne_khalil_forward
    return rne_forward(rbtdef, geom, ifunc)


def rne_backward(rbtdef, geom, fw_results, ifunc=None):
    if rbtdef._dh_convention == 'standard':
        rne_backward = rne_park_backward
    elif rbtdef._dh_convention == 'modified':
        rne_backward = rne_khalil_backward
    return rne_backward(rbtdef, geom, fw_results, ifunc)


def rne(rbtdef, geom, ifunc=None):
    '''Generate joints generic forces/torques equation.'''

    fw_results = rne_forward(rbtdef, geom, ifunc)
    invdyn = rne_backward(rbtdef, geom, fw_results, ifunc=ifunc)

    return invdyn


def gravityterm(rbtdef, geom, ifunc=None):
    '''Generate gravity force equation.'''
    if not ifunc:
        ifunc = identity
    rbtdeftmp = deepcopy(rbtdef)
    rbtdeftmp.dq = zeros((rbtdeftmp.dof, 1))
    rbtdeftmp.ddq = zeros((rbtdeftmp.dof, 1))
    geomtmp = Geometry(rbtdeftmp)
    return rne(rbtdeftmp, geomtmp, ifunc)


def coriolisterm(rbtdef, geom, ifunc=None):
    '''Generate Coriolis and centriptal forces equation.'''
    if not ifunc:
        ifunc = identity
    rbtdeftmp = deepcopy(rbtdef)
    rbtdeftmp.gravityacc = zeros((3, 1))
    rbtdeftmp.ddq = zeros((rbtdeftmp.dof, 1))
    geomtmp = Geometry(rbtdeftmp)
    return rne(rbtdeftmp, geomtmp, ifunc)


def inertiamatrix(rbtdef, geom, ifunc=None):
    '''Generate mass matrix.'''

    if not ifunc:
        ifunc = identity

    M = zeros((rbtdef.dof, rbtdef.dof))

    rbtdeftmp = deepcopy(rbtdef)
    rbtdeftmp.gravityacc = zeros((3, 1))
    rbtdeftmp.dq = zeros((rbtdeftmp.dof, 1))

    for i in range(M.rows):
        rbtdeftmp.ddq = zeros((rbtdeftmp.dof, 1))
        rbtdeftmp.ddq[i] = 1
        geomtmp = Geometry(rbtdeftmp)

        fw_results = rne_forward(rbtdeftmp, geomtmp, ifunc)
        Mcoli = rne_backward(rbtdeftmp, geomtmp, fw_results, ifunc=ifunc)

        # It's done like this since M is symmetric:
        M[:, i] = (M[i, :i].T) .col_join(Mcoli[i:, :])

    return M
