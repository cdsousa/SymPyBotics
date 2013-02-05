# -*- coding: utf-8 -*-

###############################################################################
#  SymPyBotics: Symbolic Robotics Toolbox using Python and SymPy
#
#      Copyright (C) 2012, 2013 Cristóvão Sousa <crisjss@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################


import sympy

_id = lambda x: x

class Kinematics(object):
  """Robot symbolic Jacobians."""

  def __init__(self, robot , geom, ifunc=None ):
    
    if not ifunc: ifunc = _id
    
    def sym_skew(v):
      return sympy.Matrix([[    0, -v[2],  v[1]],
                           [ v[2],     0, -v[0]],
                           [-v[1],  v[0],     0]])
    
    # extend z and p so that z[-1] and p[-1] return values from base frame
    z_ext = geom.z + [sympy.Matrix([0,0,1])]
    p_ext = geom.p + [sympy.zeros(3,1)]

    self.Jp = list(range(robot.dof))
    for l in range(robot.dof):
      self.Jp[l] = sympy.zeros((3 ,robot.dof))
      for j in range(l+1):
        if robot.links_sigma[j]:
          self.Jp[l][0:3, j] = ifunc( z_ext[j-1] )
        else:
          self.Jp[l][0:3, j] = ifunc( z_ext[j-1].cross( ( p_ext[l] - p_ext[j-1] ) ).transpose() )

    self.Jo = list(range(robot.dof))
    for l in range(robot.dof):
      self.Jo[l] = sympy.zeros((3,robot.dof))
      for j in range(l+1):
        if robot.links_sigma[j]:
          self.Jo[l][0:3, j] = sympy.zeros((3,1))
        else:
          self.Jo[l][0:3, j] = ifunc( z_ext[j-1] )

    self.J = list(range(robot.dof))
    for l in range(robot.dof):
      self.J[l] = self.Jp[l].col_join( self.Jo[l] )
      

    self.Jcp = list(range(robot.dof))
    self.Jco = self.Jo
    for l in range(robot.dof):
      self.Jcp[l] = ifunc( self.Jp[l] - sym_skew( geom.R[l]*sympy.Matrix(robot.l[l]) ) * self.Jo[l] )

    self.Jc = list(range(robot.dof))
    for l in range(robot.dof):
      self.Jc[l] = self.Jcp[l].col_join( self.Jco[l] )

  

