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


import copy
import sympy


def is_compound(expression):
  if expression.is_Atom or (-expression).is_Atom :
    return False
  else :
    return True

def intervar( expression, varpool, varrepr = '', poolrepr = '', condition_func=None, simplify_func=None ):
    if simplify_func:
      expression = simplify_func( expression )
    if condition_func and not condition_func(expression):
        return expression
    else:
        latexname = None
        if not varrepr:
            varrepr = str(len(varpool))
            if not poolrepr:
                poolrepr = 'ivar_'
        name = poolrepr+varrepr
        newvar = sympy.Symbol(poolrepr+varrepr,real=True)
        varpool.append( (newvar, expression) )
        return newvar
    
    
def m_intervar( m_exp, varpool, varrepr = '', poolrepr = '', condition_func = None, simplify_func = None ):
    m_exp_out = copy.copy(m_exp)
    repridx = None
    for e in [(i,j) for i in range(m_exp.rows) for j in range(m_exp.cols)] :
        if varrepr : repridx = varrepr+'_'+str(e[0])+'_'+str(e[1])
        m_exp_out[e] = intervar( m_exp[e] , varpool, repridx, poolrepr, condition_func, simplify_func )
    return m_exp_out


def genfunc_m_intervar( gen_intervars, ivars ):

  if gen_intervars:
    if isinstance( gen_intervars , str ):
      poolrepr = gen_intervars
    else:
      poolrepr = 'ivar_'
    def m_intervar_func( matrix, varrepr='' ):
      simplify_func = None
      #simplify_func = sympy.trigsimp
      return m_intervar( matrix, ivars, poolrepr=poolrepr, condition_func = is_compound, simplify_func = simplify_func )
  else:
    def m_intervar_func( matrix, varrepr='' ):
      return matrix

  return m_intervar_func