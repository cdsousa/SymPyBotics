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

    self.tau_ivs, self.tau = memoize(dynamic_algorithms.gen_tau_rne)( True, rbt )
    self.Y_ivs, self.Y = memoize(dynamic_algorithms.gen_regressor_rne)( True, rbt, usefricdyn=usefricdyn )
    self.M_ivs, self.M = memoize(dynamic_algorithms.gen_massmatrix_rne)( True, rbt )
    self.c_ivs, self.c = memoize(dynamic_algorithms.gen_ccfterm_rne)( True, rbt )
    self.g_ivs, self.g = memoize(dynamic_algorithms.gen_gravterm_rne)( True, rbt )
    self.f = memoize(dynamic_algorithms.gen_fricterm)( rbt )

    self.delta = sympy.Matrix( rbt.dynparms(usefricdyn=usefricdyn) )
    self.n_delta =  len( self.delta )

    self.func_def_tau = memoize(codegen_robot.dyn_matrix_to_func)( 'python', self.tau_ivs,  self.tau, 'tau_func', 2, rbt.dof, self.delta  )
    self.func_def_regressor = memoize(codegen_robot.dyn_matrix_to_func)( 'python', self.Y_ivs, self.Y, 'regressor_func', 2, rbt.dof  )
    self.func_def_M = memoize(codegen_robot.dyn_matrix_to_func)( 'python', self.M_ivs,  self.M, 'M_func', 0, rbt.dof, self.delta  )
    self.func_def_c = memoize(codegen_robot.dyn_matrix_to_func)( 'python', self.c_ivs,  self.c, 'c_func', 1, rbt.dof, self.delta  )
    self.func_def_g = memoize(codegen_robot.dyn_matrix_to_func)( 'python', self.g_ivs,  self.g, 'g_func', 0, rbt.dof, self.delta  )
    self.func_def_f = memoize(codegen_robot.dyn_matrix_to_func)( 'python', [],  self.f, 'f_func', 1, rbt.dof, self.delta  )

    global sin, cos, sign
    sin = numpy.sin
    cos = numpy.cos
    sign = numpy.sign
    compatible_exec(self.func_def_regressor,globals())
    Pb, Pd, Kd = memoize(dynamic_algorithms.find_dyn_parm_deps)( rbt.dof, self.n_delta, regressor_func )
    
    self.Pb = sympy.Matrix(Pb).applyfunc(lambda x: x.nsimplify())
    self.Pd = sympy.Matrix(Pd).applyfunc(lambda x: x.nsimplify())
    self.Kd = sympy.Matrix(Kd).applyfunc(lambda x: x.nsimplify())

    self.base_idxs = ( numpy.matrix([[i for i in range(self.n_delta)]]) * numpy.matrix(Pb) ).astype(float).astype(int).tolist()[0]
    
    self.beta = ( self.Pb.T + self.Kd * self.Pd.T ) * self.delta
    self.n_beta = len( self.beta )

