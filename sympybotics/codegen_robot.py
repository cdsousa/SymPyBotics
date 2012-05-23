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


from . import codegen
from . import geomkinem



def gen_geometric_kinematic_code( robot ):

  geomauxv, geom = geomkinem.gen_geometric_model(robot,True)
  kinemauxv, kinem = geomkinem.gen_kinematic_model(robot,geom,True)

  all_p = []
  for p in geom.pi:
    all_p += p.mat
  all_R = []
  for R in geom.Ri:
    all_R += R.mat
  all_J = []
  for i in range(len(kinem.Jpi)):
    all_J += ( kinem.Jpi[i].col_join(kinem.Joi[i]) ).mat

  geomkinem_code = codegen.optimize_code( ( geomauxv + kinemauxv, all_p + all_R + all_J ), ivarnames='aux' )

  return geomkinem_code


def _gen_q_dq_ddq_subs(dof):
  subs = []
  for i in reversed(range(dof)):
    subs.append( ( 'ddq'+str(i+1)  ,'ddq['+str(i)+']' ) )
    subs.append( ( 'dq'+str(i+1)  ,'dq['+str(i)+']' ) )
    subs.append( ( 'q'+str(i+1)  ,'q['+str(i)+']' ) )
  return subs

def _gen_parms_subs( parms_symbols, name = 'parms' ):
    subs = []
    for i in reversed(range(len(parms_symbols))):
        subs.append( ( str(parms_symbols[i]), name+'['+str(i)+']' ) )
    return subs

 
def dyn_matrix_to_func( lang, ivars, matrix, funcname, qderivlevel, dof, dynparam_symbols=[] ):
    func_parms = []
    subs_pairs = []
    if dynparam_symbols:
        func_parms.append('parms')
        subs_pairs += _gen_parms_subs(dynparam_symbols,'parms')
    if qderivlevel >= 0:
        for i in range(qderivlevel+1):
            func_parms.append('d'*i+'q')
        subs_pairs += _gen_q_dq_ddq_subs(dof)
    return codegen.sympymatrix_to_func(  lang, ivars, matrix, funcname, func_parms, subs_pairs )


