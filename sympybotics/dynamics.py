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

import copy
import sympy
import numpy

from . import geometry


def _skew(v):
  return sympy.Matrix( [ [     0, -v[2],  v[1] ],
                   [  v[2],     0, -v[0] ],
                   [ -v[1],  v[0],     0 ] ] )


def _Adj(G,g):
  R = G[0:3,0:3];  p = G[0:3,3]
  return ( R.row_join(sympy.zeros(3)) ).col_join( (_skew(p)*R).row_join(R) ) * g

def _Adjdual(G,g):
  R = G[0:3,0:3];  p = G[0:3,3]
  return ( ( R.row_join(sympy.zeros(3)) ).col_join( (_skew(p)*R).row_join(R) ) ).transpose() * g

def _adj(g,h):
  wg = g[0:3,0];  vg = g[3:6,0]
  return ( _skew(wg).row_join(sympy.zeros(3)) ).col_join( (_skew(vg)).row_join(_skew(wg)) ) * h

def _adjdual(g,h):
  wg = g[0:3,0];  vg = g[3:6,0]
  return ( ( _skew(wg).row_join(sympy.zeros(3)) ).col_join( (_skew(vg)).row_join(_skew(wg)) ) ).transpose() * h



_id = lambda x: x


def _forward_rne(rbt, ifunc=None):
  '''RNE forward pass.'''
  
  if not ifunc: ifunc = _id
  
  V = list(range(0,rbt.dof+1))
  dV = list(range(0,rbt.dof+1))

  V[-1] = sympy.zeros((6,1))
  dV[-1] = - sympy.zeros((3,1)).col_join( rbt.gravity )

  # Forward
  for i in range(rbt.dof):

    V[i] =  _Adj( rbt.geom.Tdh_inv[i], V[i-1] )  +  rbt.geom.S[i] * rbt.dq[i]
    V[i] = ifunc( V[i] )

    dV[i] =  rbt.geom.S[i] * rbt.ddq[i]  +  _Adj( rbt.geom.Tdh_inv[i], dV[i-1] )  +  _adj(  _Adj( rbt.geom.Tdh_inv[i], V[i-1] ),  rbt.geom.S[i] * rbt.dq[i] )
    dV[i] = ifunc( dV[i] )

  return V, dV


def _backward_rne(rbt, Llm, V, dV, ifunc=None):
  '''RNE backward pass.'''
  
  if not ifunc: ifunc = _id

  # extend Tdh_inv so that Tdh_inv[dof] return identity
  Tdh_inv = rbt.geom.Tdh_inv + [sympy.eye(4)]

  F = list(range(rbt.dof+1))
  F[rbt.dof] = sympy.zeros((6,1))
  
  tau = sympy.zeros((rbt.dof,1))
  
  # Backward
  for i in range(rbt.dof-1,-1,-1):

    F[i] =  _Adjdual( Tdh_inv[i+1], F[i+1] )  +  Llm[i] * dV[i]  -  _adjdual( V[i],  Llm[i] * V[i] )

    F[i] = ifunc(F[i])

    tau[i] =  ifunc(( rbt.geom.S[i].transpose() *  F[i] )[0])

  return tau




def gen_tau_rne(rbt, ifunc=None):
  '''Generate joints generic forces/torques equation.'''
  
  if not ifunc: ifunc = _id
      
  Llm = list(range(rbt.dof))

  for i in range( rbt.dof ):
    Llm[i] = (rbt.L[i].row_join(_skew(rbt.l[i])) ).col_join( (-_skew( rbt.l[i]) ).row_join(sympy.eye(3)*rbt.m[i]))

  V, dV = _forward_rne( rbt, ifunc )
  tau = _backward_rne( rbt, Llm, V, dV, ifunc )

  return tau



def gen_regressor_rne(rbt, usefricdyn=False, ifunc=None):
  '''Generate regression matrix.'''
  
  if not ifunc: ifunc = _id

  V, dV = _forward_rne( rbt, ifunc )

  dynparms = rbt.dynparms( usefricdyn = usefricdyn )

  Y = sympy.zeros( ( rbt.dof, len(dynparms) ) )
  
  if usefricdyn:
    fric = gen_fricterm(rbt)
    fric_dict = dict( zip( rbt.fc, [0]*len(rbt.fc) ) )
    fric_dict.update( dict( zip( rbt.fv, [0]*len(rbt.fv) ) ) )

  for p,parm in enumerate(dynparms):

    Llm = list(range(rbt.dof))

    for i in range(rbt.dof):
      L = rbt.L[i].applyfunc(lambda x: 1 if x == parm else 0 )
      r = rbt.l[i].applyfunc(lambda x: 1 if x == parm else 0 )
      m = 1 if rbt.m[i] == parm else 0
      
      Llm[i] = ( ( L.row_join(_skew(r)) ).col_join( (-_skew(r) ).row_join(sympy.eye(3)*m) ) )

    Y[:,p] = _backward_rne( rbt, Llm, V, dV, ifunc )

    if usefricdyn:
      select = copy.copy(fric_dict)
      select.update( { parm: 1 } )
      Y[:,p] = ifunc( Y[:,p] + fric.subs(select) )

  return Y



def gen_gravterm_rne(rbt, ifunc=None):
  '''Generate gravity force equation.'''
  if not ifunc: ifunc = _id
  rbttmp = copy.deepcopy(rbt)
  rbttmp.dq = sympy.zeros((rbttmp.dof,1))
  rbttmp.ddq = sympy.zeros((rbttmp.dof,1))
  rbttmp.geom = geometry.Geometry(rbttmp)
  return gen_tau_rne(rbttmp, ifunc)


def gen_ccfterm_rne(rbt, ifunc=None):
  '''Generate Coriolis and centriptal forces equation.'''
  if not ifunc: ifunc = _id
  rbttmp = copy.deepcopy(rbt)
  rbttmp.gravity = sympy.zeros((3,1))
  rbttmp.ddq = sympy.zeros((rbttmp.dof,1))
  rbttmp.geom = geometry.Geometry(rbttmp)
  return gen_tau_rne(rbttmp, ifunc)



def gen_massmatrix_rne(rbt, ifunc=None):
  '''Generate mass matrix.'''
  
  if not ifunc: ifunc = _id

  Llm = list(range(rbt.dof))

  for i in range( rbt.dof ):
    Llm[i] = ( rbt.L[i].row_join(_skew(rbt.l[i])) ).col_join( (-_skew( rbt.l[i]) ).row_join(sympy.eye(3)*rbt.m[i]) )

  M = sympy.zeros((rbt.dof,rbt.dof))

  rbttmp = copy.deepcopy( rbt )
  rbttmp.gravity = sympy.zeros((3,1))
  rbttmp.dq = sympy.zeros((rbttmp.dof,1))

  for i in range( M.rows ):
    rbttmp.ddq = sympy.zeros((rbttmp.dof,1))
    rbttmp.ddq[i] = 1
    rbttmp.geom = geometry.Geometry(rbttmp)

    V, dV = _forward_rne( rbttmp, ifunc )
    Mcoli = _backward_rne( rbttmp, Llm, V, dV, ifunc )

    # It's done like this since M is symmetric:
    M[:,i] = ( M[i,:i].T ) .col_join( Mcoli[i:,:] )

  return M



def gen_fricterm(rbt, ifunc=None):
  '''Generate friction forces (simple Coulomb and viscouse model).'''
  if not ifunc: ifunc = _id
  fric = sympy.zeros((rbt.dof,1))
  for i in range( rbt.dof ):
    fric[i] = ifunc( rbt.fv[i] * rbt.dq[i] + rbt.fc[i] * sympy.sign(rbt.dq[i]) )
  return fric



def find_dyn_parm_deps( dof, parm_num, regressor_func ):
  '''Find dynamic parameter dependencies (i.e., regressor column dependencies).'''

  samples=10000
  round = 10

  pi = numpy.pi

  Z = numpy.zeros( ( dof*samples, parm_num ) )

  for i in range(samples):
    q = [ float( numpy.random.random()*2.0*pi - pi ) for j in range(dof) ]
    dq = [ float( numpy.random.random()*2.0*pi - pi ) for j in range(dof) ]
    ddq = [ float( numpy.random.random()*2.0*pi - pi ) for j in range(dof) ]
    Z[ i*dof : i*dof+dof , : ] = numpy.matrix( regressor_func( q, dq, ddq ) ).reshape( dof, parm_num )

  R1_diag=  numpy.linalg.qr( Z, mode='economic' ).diagonal().round(round)
  dbi = []
  ddi = []
  for i,e in enumerate(R1_diag):
      if e != 0:
          dbi.append(i)
      else:
          ddi.append(i)
  dbn = len(dbi)

  P = numpy.mat(numpy.eye(parm_num))[ :, dbi+ddi ]
  Pb = P[:,:dbn]
  Pd = P[:,dbn:]

  Rbd1 = numpy.mat(numpy.linalg.qr( Z*P, mode='r' ))
  Rb1 = Rbd1[:dbn,:dbn]
  Rd1 = Rbd1[:dbn,dbn:]

  Kd = numpy.mat( ( numpy.linalg.inv(Rb1) * Rd1 ).round(round) )

  return Pb,Pd,Kd



class Dynamics(object):
  """Robot dynamic model in code form."""

  def __init__(self, rbt, usefricdyn, ifunc=None):
    
    self.usefricdyn = usefricdyn

    self.delta = sympy.Matrix( rbt.dynparms(usefricdyn=self.usefricdyn) )
    self.n_delta =  len( self.delta )

    self.tau = gen_tau_rne(rbt, ifunc)
    self.regressor = gen_regressor_rne(rbt, usefricdyn=self.usefricdyn, ifunc=ifunc)
    self.M = gen_massmatrix_rne(rbt, ifunc)
    self.c = gen_ccfterm_rne(rbt, ifunc)
    self.g = gen_gravterm_rne(rbt, ifunc)
    if self.usefricdyn:
      self.f = gen_fricterm(rbt, ifunc)
      
      
  #def find_base_parameters(self, rbt):
      
    #func_def_regressor = codegen_robot.dyn_code_to_func( 'python', self.tau.regressor_code, 'regressor_func', 2, rbt.dof  )
    #global sin, cos, sign
    #sin = numpy.sin
    #cos = numpy.cos
    #sign = numpy.sign
    #compatible_exec(func_def_regressor,globals())
    
    #Pb, Pd, Kd = find_dyn_parm_deps( rbt.dof, self.n_delta, regressor_func )
    
    #self.Pb = sympy.Matrix(Pb).applyfunc(lambda x: x.nsimplify())
    #self.Pd = sympy.Matrix(Pd).applyfunc(lambda x: x.nsimplify())
    #self.Kd = sympy.Matrix(Kd).applyfunc(lambda x: x.nsimplify())

    #self.base_idxs = ( numpy.matrix([[i for i in range(self.n_delta)]]) * numpy.matrix(Pb) ).astype(float).astype(int).tolist()[0]
    
    #self.beta = ( self.Pb.T + self.Kd * self.Pd.T ) * self.delta
    #self.n_beta = len( self.beta )
    
    #self.base_regressor = self.regressor * self.Pb
