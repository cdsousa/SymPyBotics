
import sys
import sympy
import numpy

import symcode

from . import geometry
from . import kinematics
from . import dynamics

def _fprint(x):
  print(x)
  sys.stdout.flush()

class RobotAllSymb(object):
  """Robot geometric, kinematic, and dynamic models in single symbolic expressions."""

  def __init__(self, rbtdef):
    
    self.rbtdef = rbtdef
    self.dof = rbtdef.dof
    
    self.geom = geometry.Geometry(self.rbtdef)
    self.kin = kinematics.Kinematics(self.rbtdef, self.geom)
    
    self.dyn = dynamics.Dynamics(self.rbtdef, self.geom)
    self.dyn.gen_all()

  
class RobotDynCode(object):
  """Robot dynamic model in code form."""
  
  def __init__(self, rbtdef):
    
    self.rbtdef = rbtdef
    self.dof = rbtdef.dof
    
    _fprint('generating geometric model')
    self.geom = geometry.Geometry(self.rbtdef)
    
    _fprint('generating kinematic model')
    self.kin = kinematics.Kinematics(self.rbtdef, self.geom)
    
    self.dyn = dynamics.Dynamics(self.rbtdef, self.geom)
    
    _fprint('generating tau code')
    tau_se = symcode.subexprs.Subexprs()
    self.dyn.gen_tau(tau_se.collect)
    self.tau_code  = tau_se.get(self.dyn.tau)
    
    _fprint('generating gravity term code')
    g_se = symcode.subexprs.Subexprs()
    self.dyn.gen_gravityterm(g_se.collect)
    self.g_code  = g_se.get(self.dyn.gravityterm)
    
    _fprint('generating coriolis term code')
    c_se = symcode.subexprs.Subexprs()
    self.dyn.gen_coriolisterm(c_se.collect)
    self.c_code  = c_se.get(self.dyn.coriolisterm)
    
    _fprint('generating inertia matrix code')
    M_se = symcode.subexprs.Subexprs()
    self.dyn.gen_inertiamatrix(M_se.collect)
    self.M_code  = M_se.get(self.dyn.inertiamatrix)
    
    _fprint('generating regressor matrix code')
    H_se = symcode.subexprs.Subexprs()
    self.dyn.gen_regressor(H_se.collect)
    self.H_code  = H_se.get(self.dyn.regressor)
    self._H_se =H_se
    
    self._codes = ['tau_code', 'g_code', 'c_code', 'M_code', 'H_code']
    
    if self.rbtdef.frictionmodel != None:    
      _fprint('generating friction term code')
      f_se = symcode.subexprs.Subexprs()
      self.dyn.gen_frictionterm(f_se.collect)
      self.f_code  = f_se.get(self.dyn.frictionterm)
      self._codes.append('f_code')
      
    _fprint('done')
    
    
  def calc_base_parms(self):
    
    q_subs = {q:'q[%d]'%i for i, q in enumerate(self.rbtdef.q)}
    q_subs.update({dq:'dq[%d]'%i for i, dq in enumerate(self.rbtdef.dq)})
    q_subs.update({ddq:'ddq[%d]'%i for i, ddq in enumerate(self.rbtdef.ddq)})
    func_def_regressor = symcode.generation.code_to_func('python', self.H_code, 'regressor_func', ['q','dq','ddq'], q_subs)
    global sin, cos, sign
    sin = numpy.sin
    cos = numpy.cos
    sign = numpy.sign
    exec(func_def_regressor)
    
    _fprint('calculating base parameters and regressor code')
    self.dyn.calc_base_parms(regressor_func)
    self.Hb_code  = self._H_se.get(self.H_code[1]*self.Pb)
    
    self._codes.append('Hb_code')
    
    _fprint('done')