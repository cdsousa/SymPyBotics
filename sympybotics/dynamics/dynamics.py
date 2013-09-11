import sympy
import numpy

from ..geometry import Geometry

from .rne import rne, gravityterm, coriolisterm, inertiamatrix
from .regressor import regressor
from .simple_fric import frictionterm
from .dyn_parm_dep import find_dyn_parm_deps


class Dynamics(object):

    def __init__(self, rbtdef, geom):

        self.rbtdef = rbtdef
        self.geom = geom
        self.dof = rbtdef.dof

        self.dynparms = sympy.Matrix(rbtdef.dynparms())
        self.n_dynparms = len(self.dynparms)

    def gen_invdyn(self, ifunc=None):
        self.invdyn = rne(self.rbtdef, self.geom, ifunc)

    def gen_gravityterm(self, ifunc=None):
        self.g = gravityterm(self.rbtdef, self.geom, ifunc)

    def gen_coriolisterm(self, ifunc=None):
        self.c = coriolisterm(
            self.rbtdef, self.geom, ifunc)

    def gen_inertiamatrix(self, ifunc=None):
        self.M = inertiamatrix(
            self.rbtdef, self.geom, ifunc)

    def gen_regressor(self, ifunc=None):
        self.H = regressor(self.rbtdef, self.geom, ifunc)

    def gen_frictionterm(self, ifunc=None):
        if self.rbtdef.frictionmodel == 'simple':
            self.f = frictionterm(self.rbtdef, ifunc)
        else:
            self.f = sympy.zeros(self.dof, 1)

    def gen_all(self, ifunc=None):
        self.gen_invdyn(ifunc)
        self.gen_gravityterm(ifunc)
        self.gen_coriolisterm(ifunc)
        self.gen_inertiamatrix(ifunc)
        self.gen_regressor(ifunc)

    def calc_base_parms(self, regressor_func):

        Pb, Pd, Kd = find_dyn_parm_deps(
            self.dof, self.n_dynparms, regressor_func)

        self.Pb = sympy.Matrix(Pb).applyfunc(lambda x: x.nsimplify())
        self.Pd = sympy.Matrix(Pd).applyfunc(lambda x: x.nsimplify())
        self.Kd = sympy.Matrix(Kd).applyfunc(lambda x: x.nsimplify())

        self.base_idxs = \
            (numpy.matrix([[i for i in range(self.n_dynparms)]]) *
             numpy.matrix(Pb)).astype(float).astype(int).tolist()[0]

        self.baseparms = (self.Pb.T + self.Kd * self.Pd.T) * self.dynparms
        self.n_base = len(self.baseparms)
