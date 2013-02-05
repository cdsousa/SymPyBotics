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

from . import geometry
from . import kinematics
from . import dynamics
from . import code

class RobotAllSymb(object):
  """Robot geometric, kinematic, and dynamic models in single symbolic expressions."""

  def __init__(self, rbtdef, usefricdyn=False):
    
    self.rbtdef = rbtdef
    self.dof = rbtdef.dof
    self.usefricdyn = usefricdyn
    
    self.geom = geometry.Geometry(self.rbtdef)
    self.kin = geometry.Geometry(self.rbtdef, self.geom)
    
    self.dyn = dynamics.Dynamics(self.rbtdef, self.geom, self.usefricdyn)
    self.dyn.gen_all()

  
class RobotDynCode(object):
  """Robot dynamic model in code form."""
  
  def __init__(self, rbtdef, usefricdyn=False):
    
    self.rbtdef = rbtdef
    self.dof = rbtdef.dof
    self.usefricdyn = usefricdyn
    
    self.geom = geometry.Geometry(self.rbtdef)
    self.kin = kinematics.Kinematics(self.rbtdef, self.geom)
    
    self.dyn = dynamics.Dynamics(self.rbtdef, self.geom, self.usefricdyn)
    
    tau_se = code.subexprs.Subexprs()
    g_se = code.subexprs.Subexprs()
    c_se = code.subexprs.Subexprs()
    M_se = code.subexprs.Subexprs()
    H_se = code.subexprs.Subexprs()
    
    self.dyn.gen_tau(tau_se.collect)
    self.dyn.gen_gravterm(g_se.collect)
    self.dyn.gen_ccfterm(c_se.collect)
    self.dyn.gen_massmatrix(M_se.collect)
    self.dyn.gen_regressor(H_se.collect)
    self.dyn.gen_base_parms()
        
    self.tau_code  = (tau_se.subexprs, self.dyn.tau)
    self.g_code  = (g_se.subexprs, self.dyn.g)
    self.c_code  = (c_se.subexprs, self.dyn.c)
    self.M_code  = (M_se.subexprs, self.dyn.M)
    self.H_code  = (H_se.subexprs, self.dyn.H)
    self.Hb_code  = (H_se.subexprs, self.dyn.Hb)
    
    