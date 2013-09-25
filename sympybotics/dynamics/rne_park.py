from sympy import zeros, eye
from ..utils import sym_skew as skew
from ..utils import identity


def Adj(G, g):
    R = G[0:3, 0:3]
    p = G[0:3, 3]
    return (R.row_join(zeros(3))).col_join(
        (skew(p) * R).row_join(R)) * g


def Adjdual(G, g):
    R = G[0:3, 0:3]
    p = G[0:3, 3]
    return ((R.row_join(zeros(3))).col_join(
        (skew(p) * R).row_join(R))).transpose() * g


def adj(g, h):
    wg = g[0:3, 0]
    vg = g[3:6, 0]
    return (skew(wg).row_join(zeros(3))).col_join(
        (skew(vg)).row_join(skew(wg))) * h


def adjdual(g, h):
    wg = g[0:3, 0]
    vg = g[3:6, 0]
    return ((skew(wg).row_join(zeros(3))).col_join(
        (skew(vg)).row_join(skew(wg)))).transpose() * h


def rne_park_forward(rbtdef, geom, ifunc=None):
    '''RNE forward pass.'''

    if not ifunc:
        ifunc = identity

    V = list(range(0, rbtdef.dof + 1))
    dV = list(range(0, rbtdef.dof + 1))

    V[-1] = zeros((6, 1))
    dV[-1] = - zeros((3, 1)).col_join(rbtdef.gravityacc)

    # Forward
    for i in range(rbtdef.dof):

        V[i] = ifunc(Adj(geom.Tdh_inv[i], V[i - 1])) + \
            ifunc(geom.S[i] * rbtdef.dq[i])
        V[i] = ifunc(V[i])

        dV[i] = ifunc(geom.S[i] * rbtdef.ddq[i]) + \
            ifunc(Adj(geom.Tdh_inv[i], dV[i - 1])) + \
            ifunc(adj(ifunc(Adj(geom.Tdh_inv[i], V[i - 1])),
                  ifunc(geom.S[i] * rbtdef.dq[i])))
        dV[i] = ifunc(dV[i])

    return V, dV


def rne_park_backward(rbtdef, geom, fw_results, ifunc=None):
    '''RNE backward pass.'''

    V, dV = fw_results

    if not ifunc:
        ifunc = identity

    # extend Tdh_inv so that Tdh_inv[dof] return identity
    Tdh_inv = geom.Tdh_inv + [eye(4)]

    F = list(range(rbtdef.dof + 1))
    F[rbtdef.dof] = zeros((6, 1))

    tau = zeros((rbtdef.dof, 1))

    # Backward
    for i in range(rbtdef.dof - 1, -1, -1):

        Llm = (rbtdef.L[i].row_join(skew(rbtdef.l[i]))).col_join(
            (-skew(rbtdef.l[i])).row_join(eye(3) * rbtdef.m[i]))

        F[i] = Adjdual(Tdh_inv[i + 1], F[i + 1]) + \
            Llm * dV[i] - adjdual(V[i],  Llm * V[i])
        F[i] = ifunc(F[i])

        tau[i] = ifunc((geom.S[i].transpose() * F[i])[0])

    return tau

