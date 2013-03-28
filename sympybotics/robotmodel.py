
import sys
import sympy
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
  
  def __init__(self, rbtdef, codecollectmode='unique_ops'):
    
    self.rbtdef = rbtdef
    self.dof = rbtdef.dof
    
    self.geom = geometry.Geometry(self.rbtdef)
    self.kin = kinematics.Kinematics(self.rbtdef, self.geom)
    
    self.dyn = dynamics.Dynamics(self.rbtdef, self.geom)
    
    _fprint('generating tau code')
    tau_se = symcode.subexprs.Subexprs(codecollectmode)
    self.dyn.gen_tau(tau_se.collect)
    self.tau_code  = (tau_se.subexprs, self.dyn.tau)
    
    _fprint('generating gravity term code')
    g_se = symcode.subexprs.Subexprs(codecollectmode)
    self.dyn.gen_gravityterm(g_se.collect)
    self.g_code  = (g_se.subexprs, self.dyn.gravityterm)
    
    _fprint('generating coriolis term code')
    c_se = symcode.subexprs.Subexprs(codecollectmode)
    self.dyn.gen_coriolisterm(c_se.collect)
    self.c_code  = (c_se.subexprs, self.dyn.coriolisterm)
    
    _fprint('generating inertia matrix code')
    M_se = symcode.subexprs.Subexprs(codecollectmode)
    self.dyn.gen_inertiamatrix(M_se.collect)
    self.M_code  = (M_se.subexprs, self.dyn.inertiamatrix)
    
    _fprint('generating regressor matrix code')
    H_se = symcode.subexprs.Subexprs(codecollectmode)
    self.dyn.gen_regressor(H_se.collect)
    self.H_code  = (H_se.subexprs, self.dyn.regressor)
    
    self._codes = ['tau_code', 'g_code', 'c_code', 'M_code', 'H_code']
    
    if self.rbtdef.frictionmodel != None:    
      _fprint('generating friction term code')
      f_se = symcode.subexprs.Subexprs(codecollectmode)
      self.dyn.gen_frictionterm(f_se.collect)
      self.f_code  = (f_se.subexprs, self.dyn.frictionterm)
      self._codes.append('f_code')
      
    _fprint('done')
    
    
  def calc_base_parms(self):
    
    _fprint('calculating base parameters and regressor code')
    self.dyn.calc_base_parms()
    self.Hb_code  = (self.H_code[0], self.dyn.base_regressor)
    
    self._codes.append('Hb_code')
    
    _fprint('done')


  def optimize_code(self, mode='light'):
    
    if mode == 'light':
      def optimize( code ):
        code = symcode.optimization.optim_dce_sup(code)
        code = symcode.optimization.optim_cp(code)
        return code
    elif mode == 'heavy':
      def optimize( code ):
        code = symcode.optimization.fully_optimize_code(code)
        return code
    else:
      raise Exception("No '%s' optimization mode known."%mode)
      
    for codename in self._codes:
      _fprint('optimizing %s'%codename)
      oldcode = getattr(self,codename)
      newcode = optimize(oldcode)
      setattr(self,codename,newcode)
    
    _fprint('done')
    