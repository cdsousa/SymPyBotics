from sympy import zeros, eye, Matrix
from .extra_dyn import frictionforce, driveinertiaterm
from ..utils import sym_skew as skew
from ..utils import identity


def rne_khalil_forward(rbtdef, geom, ifunc=None):
    '''RNE forward pass.'''

    if not ifunc:
        ifunc = identity

    w = list(range(0, rbtdef.dof + 1))
    dw = list(range(0, rbtdef.dof + 1))
    dV = list(range(0, rbtdef.dof + 1))
    U = list(range(0, rbtdef.dof + 1))

    w[-1] = zeros((3, 1))
    dw[-1] = zeros((3, 1))
    dV[-1] = -rbtdef.gravityacc
    U[-1] = zeros((3, 3))

    z = Matrix([0, 0, 1])

    # Forward
    for i in range(rbtdef.dof):

        s = rbtdef._links_sigma[i]
        ns = 1 - s

        w_pj = geom.Rdh[i].T * w[i - 1]

        w[i] = w_pj + ns * rbtdef.dq[i] * z
        w[i] = ifunc(w[i])

        dw[i] = geom.Rdh[i].T * dw[i - 1] + ns * \
            (rbtdef.ddq[i] * z + w_pj.cross(rbtdef.dq[i] * z).T)
        dw[i] = ifunc(dw[i])

        dV[i] = geom.Rdh[i].T * (dV[i - 1] + U[i - 1] * geom.pdh[i]) + s * (
            rbtdef.ddq[i] * z + 2 * w_pj.cross(rbtdef.dq[i] * z).T)
        dV[i] = ifunc(dV[i])

        U[i] = skew(dw[i]) + skew(w[i]) ** 2
        U[i] = ifunc(U[i])

    return w, dw, dV, U


def rne_khalil_backward(rbtdef, geom, fw_results, ifunc=None):
    '''RNE backward pass.'''

    w, dw, dV, U = fw_results

    if not ifunc:
        ifunc = identity

    # extend Rdh so that Rdh[dof] return identity
    Rdh = geom.Rdh + [eye(3)]
    # extend pdh so that pRdh[dof] return zero
    pdh = geom.pdh + [zeros((3, 1))]

    F = list(range(rbtdef.dof))
    M = list(range(rbtdef.dof))
    f = list(range(rbtdef.dof + 1))
    m = list(range(rbtdef.dof + 1))

    f[rbtdef.dof] = zeros((3, 1))
    m[rbtdef.dof] = zeros((3, 1))

    z = Matrix([0, 0, 1])

    tau = zeros((rbtdef.dof, 1))

    fric = frictionforce(rbtdef)
    Idrive = driveinertiaterm(rbtdef)

    # Backward
    for i in range(rbtdef.dof - 1, -1, -1):

        s = rbtdef._links_sigma[i]
        ns = 1 - s

        F[i] = rbtdef.m[i] * dV[i] + U[i] * Matrix(rbtdef.l[i])
        F[i] = ifunc(F[i])

        M[i] = rbtdef.L[i] * dw[i] + w[i].cross(
            rbtdef.L[i] * w[i]).T + Matrix(rbtdef.l[i]).cross(dV[i]).T
        M[i] = ifunc(M[i])

        f_nj = Rdh[i + 1] * f[i + 1]

        f[i] = F[i] + f_nj  # + f_e[i]
        f[i] = ifunc(f[i])

        m[i] = M[i] + Rdh[i + 1] * m[i + 1] + \
            pdh[i + 1].cross(f_nj).T  # + m_e[i]
        m[i] = ifunc(m[i])

        tau[i] = ifunc(((s * f[i] + ns * m[i]).T * z)[0] + fric[i] + Idrive[i])

    return tau
