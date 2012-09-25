# -*- coding: utf-8 -*-

###############################################################################
#  SymPyBotics: Symbolic Robotics Toolbox using Python and SymPy
#
#      Copyright (C) 2012 Cristóvão Sousa <crisjss@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

##### For both Python 2 and Python 3 compatiblility:
import sys
if sys.version_info.major >= 3:
    def compatible_exec(source, g=None,l=None):
      if g:
        if l:
          exec(source, g, l)
        else:
          exec(source, g)
      else:
        exec(source)
else:
    eval(compile("""\
def compatible_exec(source, g=None,l=None):
      if g:
        if l:
          exec source in g, l
        else:
          exec source in g
      else:
        exec source
""",
    "<exec_function>", "exec"))

import numpy
import sympy
from . import dynamic_algorithms
from . import codegen
from . import codegen_robot
from . import tools


class Dyn(object):
  """Robot dynamic model in symbolic and python code forms."""

  def __init__( self, rbt, usefricdyn, memoize_func=None ):

    if memoize_func:
      memoize = memoize_func
    else:
      def memoize( func, extra_deps='' ):
        def decorated_function( *args, **kwargs ):
          return func( *args, **kwargs )
        return decorated_function

    self.dof = rbt.dof

    self.delta = sympy.Matrix( rbt.dynparms(usefricdyn=usefricdyn) )
    self.n_delta =  len( self.delta )

    tau_ivs, tau = memoize(dynamic_algorithms.gen_tau_rne)( True, rbt )
    regressor_ivs, regressor = memoize(dynamic_algorithms.gen_regressor_rne)( True, rbt, usefricdyn=usefricdyn )
    M_ivs, M = memoize(dynamic_algorithms.gen_massmatrix_rne)( True, rbt )
    c_ivs, c = memoize(dynamic_algorithms.gen_ccfterm_rne)( True, rbt )
    g_ivs, g = memoize(dynamic_algorithms.gen_gravterm_rne)( True, rbt )
    if usefricdyn:
      f = memoize(dynamic_algorithms.gen_fricterm)( rbt )

    self.tau_code = memoize(codegen.optimize_code)( (tau_ivs, tau.mat), ivarnames='aux' )
    self.regressor_code = memoize(codegen.optimize_code)( (regressor_ivs, regressor.mat), ivarnames='aux' )
    self.M_code = memoize(codegen.optimize_code)( (M_ivs, M.mat), ivarnames='aux' )
    self.c_code = memoize(codegen.optimize_code)( (c_ivs, c.mat), ivarnames='aux' )
    self.g_code = memoize(codegen.optimize_code)( (g_ivs, g.mat), ivarnames='aux' )
    if usefricdyn:
      self.f_code = memoize(codegen.optimize_code)( ([], f.mat), ivarnames='aux' )
    
    func_def_regressor = memoize(codegen_robot.dyn_code_to_func)( 'python', self.regressor_code, 'regressor_func', 2, rbt.dof  )
    global sin, cos, sign
    sin = numpy.sin
    cos = numpy.cos
    sign = numpy.sign
    compatible_exec(func_def_regressor,globals())
    Pb, Pd, Kd = memoize(dynamic_algorithms.find_dyn_parm_deps)( rbt.dof, self.n_delta, regressor_func )
    
    self.Pb = sympy.Matrix(Pb).applyfunc(lambda x: x.nsimplify())
    self.Pd = sympy.Matrix(Pd).applyfunc(lambda x: x.nsimplify())
    self.Kd = sympy.Matrix(Kd).applyfunc(lambda x: x.nsimplify())

    self.base_idxs = ( numpy.matrix([[i for i in range(self.n_delta)]]) * numpy.matrix(Pb) ).astype(float).astype(int).tolist()[0]
    
    self.beta = ( self.Pb.T + self.Kd * self.Pd.T ) * self.delta
    self.n_beta = len( self.beta )

    self.gen_member_funcs()


  def gen_member_funcs(self):
    
    tau_str = codegen_robot.dyn_code_to_func( 'python', self.tau_code, 'tau', 2, self.dof, self.delta )
    regressor_str = codegen_robot.dyn_code_to_func( 'python', self.regressor_code, 'regressor', 2, self.dof, self.delta )
    M_str = codegen_robot.dyn_code_to_func( 'python', self.M_code, 'M', 0, self.dof )
    g_str = codegen_robot.dyn_code_to_func( 'python', self.c_code, 'c', 1, self.dof, self.delta )
    c_str = codegen_robot.dyn_code_to_func( 'python', self.g_code, 'g', 0, self.dof, self.delta )
    if hasattr(self, 'f_code'):
      f_str = codegen_robot.dyn_code_to_func( 'python', self.f_code, 'f', 1, self.dof, self.delta )

    funcs = [ ('tau', self.dof, 1),
              ('regressor', self.dof, self.n_delta),
              ('M', self.dof, self.dof),
              ('c', self.dof, 1),
              ('g', self.dof, 1)
            ]
            
    if hasattr(self, 'f_code'):
      funcs.append( ('f', self.dof, 1 ) )
    
    for func,m,n in funcs:
      func_str = eval(func+'_str')
      func_lines = func_str.splitlines()
      func_return = func_lines[-1].split()[1]
      func_lines[-1] = func_lines[-1].replace( func_return, 'sympy.Matrix( ' + str(m) + ', ' + str(n) + ', ' + func_return + ' )' )
      func_str = '\n'.join(func_lines)
      compatible_exec( func_str, globals() )

    self.tau = tau
    self.regressor = regressor
    self.M = M
    self.c = c
    self.g = g
    if hasattr(self, 'f_code'):
      self.f = f

    