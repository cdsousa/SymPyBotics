""" Utilities """

from sympy import Matrix


identity = lambda x: x


def sym_skew(v):
    return Matrix([[0,     -v[2],  v[1]],
                   [v[2],      0, -v[0]],
                   [-v[1],  v[0],     0]])
