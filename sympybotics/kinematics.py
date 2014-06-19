
import sympy

_id = lambda x: x


class Kinematics(object):

    """Robot symbolic Jacobians.

        kinobj.J: list of link frame Jacobians - complete (6 x N):
                   [linear_velocity
                    angular_velocity] = J * joint_velocities
        kinobj.Jc: list of link center-of-mass Jacobians - complete
        kinobj.Jp: list of link frame Jacobians - linear velocity part only
        kinobj.Jo: list of link frame Jacobians - angular velocity part only
        kinobj.Jcp: list of link center-of-mass Jacobians - linear part
        kinobj.Jco: list of link center-of-mass Jacobians - angular part
    """

    def __init__(self, robotdef, geom, ifunc=None):

        if not ifunc:
            ifunc = _id

        self.rbtdef = robotdef
        self.geom = geom
        self.dof = self.rbtdef.dof

        def sym_skew(v):
            return sympy.Matrix([[0,    -v[2],  v[1]],
                                 [v[2],     0, -v[0]],
                                 [-v[1], v[0],     0]])

        if self.rbtdef._dh_convention == 'standard':

            # extend z and p so that z[-1] and p[-1] return values from base
            # frame
            z_ext = geom.z + [sympy.Matrix([0, 0, 1])]
            p_ext = geom.p + [sympy.zeros(3, 1)]

            self.Jp = list(range(self.rbtdef.dof))
            for l in range(self.rbtdef.dof):
                self.Jp[l] = sympy.zeros(3, self.rbtdef.dof)
                for j in range(l + 1):
                    if self.rbtdef._links_sigma[j]:
                        self.Jp[l][0:3, j] = ifunc(z_ext[j - 1])
                    else:
                        self.Jp[l][0:3, j] = ifunc(z_ext[j - 1].cross(
                            (p_ext[l] - p_ext[j - 1])).reshape(3, 1))

            self.Jo = list(range(self.rbtdef.dof))
            for l in range(self.rbtdef.dof):
                self.Jo[l] = sympy.zeros(3, self.rbtdef.dof)
                for j in range(l + 1):
                    if self.rbtdef._links_sigma[j]:
                        self.Jo[l][0:3, j] = sympy.zeros(3, 1)
                    else:
                        self.Jo[l][0:3, j] = ifunc(z_ext[j - 1])

        elif self.rbtdef._dh_convention == 'modified':

            self.Jp = list(range(self.rbtdef.dof))
            for l in range(self.rbtdef.dof):
                self.Jp[l] = sympy.zeros(3, self.rbtdef.dof)
                for j in range(l + 1):
                    if self.rbtdef._links_sigma[j]:
                        self.Jp[l][0:3, j] = ifunc(geom.z[j])
                    else:
                        self.Jp[l][0:3, j] = ifunc(geom.z[j].cross(
                            (geom.p[l] - geom.p[j])).reshape(3, 1))

            self.Jo = list(range(self.rbtdef.dof))
            for l in range(self.rbtdef.dof):
                self.Jo[l] = sympy.zeros(3, self.rbtdef.dof)
                for j in range(l + 1):
                    if self.rbtdef._links_sigma[j]:
                        self.Jo[l][0:3, j] = sympy.zeros(3, 1)
                    else:
                        self.Jo[l][0:3, j] = ifunc(geom.z[j])

        self.J = list(range(self.rbtdef.dof))
        for l in range(self.rbtdef.dof):
            self.J[l] = self.Jp[l].col_join(self.Jo[l])

        self.Jcp = list(range(self.rbtdef.dof))
        self.Jco = self.Jo
        for l in range(self.rbtdef.dof):
            self.Jcp[l] = ifunc(self.Jp[l] - sym_skew(
                geom.R[l] * sympy.Matrix(self.rbtdef.l[l])) * self.Jo[l])

        self.Jc = list(range(self.rbtdef.dof))
        for l in range(self.rbtdef.dof):
            self.Jc[l] = self.Jcp[l].col_join(self.Jco[l])
