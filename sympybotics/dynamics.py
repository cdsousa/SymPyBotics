
import copy
import sympy
import numpy

from . import geometry
from . import symcode


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


def _forward_rne(rbtdef, geom, ifunc=None):
  '''RNE forward pass.'''
  
  if not ifunc: ifunc = _id
  
  V = list(range(0,rbtdef.dof+1))
  dV = list(range(0,rbtdef.dof+1))

  V[-1] = sympy.zeros((6,1))
  dV[-1] = - sympy.zeros((3,1)).col_join( rbtdef.gravity )

  # Forward
  for i in range(rbtdef.dof):

    V[i] =  _Adj( geom.Tdh_inv[i], V[i-1] )  +  geom.S[i] * rbtdef.dq[i]
    V[i] = ifunc( V[i] )

    dV[i] =  geom.S[i] * rbtdef.ddq[i]  +  _Adj( geom.Tdh_inv[i], dV[i-1] )  +  _adj(  _Adj( geom.Tdh_inv[i], V[i-1] ),  geom.S[i] * rbtdef.dq[i] )
    dV[i] = ifunc( dV[i] )

  return V, dV


def _backward_rne(rbtdef, geom, Llm, V, dV, ifunc=None):
  '''RNE backward pass.'''
  
  if not ifunc: ifunc = _id

  # extend Tdh_inv so that Tdh_inv[dof] return identity
  Tdh_inv = geom.Tdh_inv + [sympy.eye(4)]

  F = list(range(rbtdef.dof+1))
  F[rbtdef.dof] = sympy.zeros((6,1))
  
  tau = sympy.zeros((rbtdef.dof,1))
  
  # Backward
  for i in range(rbtdef.dof-1,-1,-1):

    F[i] =  _Adjdual( Tdh_inv[i+1], F[i+1] )  +  Llm[i] * dV[i]  -  _adjdual( V[i],  Llm[i] * V[i] )

    F[i] = ifunc(F[i])

    tau[i] =  ifunc(( geom.S[i].transpose() *  F[i] )[0])

  return tau




def _gen_tau_rne(rbtdef, geom, ifunc=None):
  '''Generate joints generic forces/torques equation.'''
  
  if not ifunc: ifunc = _id
      
  Llm = list(range(rbtdef.dof))

  for i in range( rbtdef.dof ):
    Llm[i] = (rbtdef.L[i].row_join(_skew(rbtdef.l[i])) ).col_join( (-_skew( rbtdef.l[i]) ).row_join(sympy.eye(3)*rbtdef.m[i]))

  V, dV = _forward_rne( rbtdef, geom, ifunc )
  tau = _backward_rne( rbtdef, geom, Llm, V, dV, ifunc )

  return tau



def _gen_regressor_rne(rbtdef, geom, ifunc=None):
  '''Generate regression matrix.'''
  
  if not ifunc: ifunc = _id

  V, dV = _forward_rne( rbtdef, geom, ifunc )

  dynparms = rbtdef.dynparms()

  Y = sympy.zeros( ( rbtdef.dof, len(dynparms) ) )
  
  if rbtdef.frictionmodel == 'simple':
    fric = _gen_fricterm(rbtdef)
    fric_dict = dict( zip( rbtdef.fc, [0]*len(rbtdef.fc) ) )
    fric_dict.update( dict( zip( rbtdef.fv, [0]*len(rbtdef.fv) ) ) )

  for p,parm in enumerate(dynparms):

    Llm = list(range(rbtdef.dof))

    for i in range(rbtdef.dof):
      L = rbtdef.L[i].applyfunc(lambda x: 1 if x == parm else 0 )
      r = rbtdef.l[i].applyfunc(lambda x: 1 if x == parm else 0 )
      m = 1 if rbtdef.m[i] == parm else 0
      
      Llm[i] = ( ( L.row_join(_skew(r)) ).col_join( (-_skew(r) ).row_join(sympy.eye(3)*m) ) )

    Y[:,p] = _backward_rne( rbtdef, geom, Llm, V, dV, ifunc )

    if rbtdef.frictionmodel == 'simple':
      select = copy.copy(fric_dict)
      select.update( { parm: 1 } )
      Y[:,p] = ifunc( Y[:,p] + fric.subs(select) )

  return Y



def _gen_gravterm_rne(rbtdef, geom, ifunc=None):
  '''Generate gravity force equation.'''
  if not ifunc: ifunc = _id
  rbtdeftmp = copy.deepcopy(rbtdef)
  rbtdeftmp.dq = sympy.zeros((rbtdeftmp.dof,1))
  rbtdeftmp.ddq = sympy.zeros((rbtdeftmp.dof,1))
  geomtmp = geometry.Geometry(rbtdeftmp)
  return _gen_tau_rne(rbtdeftmp, geomtmp, ifunc)


def _gen_ccfterm_rne(rbtdef, geom, ifunc=None):
  '''Generate Coriolis and centriptal forces equation.'''
  if not ifunc: ifunc = _id
  rbtdeftmp = copy.deepcopy(rbtdef)
  rbtdeftmp.gravity = sympy.zeros((3,1))
  rbtdeftmp.ddq = sympy.zeros((rbtdeftmp.dof,1))
  geomtmp = geometry.Geometry(rbtdeftmp)
  return _gen_tau_rne(rbtdeftmp, geomtmp, ifunc)



def _gen_massmatrix_rne(rbtdef, geom, ifunc=None):
  '''Generate mass matrix.'''
  
  if not ifunc: ifunc = _id

  Llm = list(range(rbtdef.dof))

  for i in range( rbtdef.dof ):
    Llm[i] = ( rbtdef.L[i].row_join(_skew(rbtdef.l[i])) ).col_join( (-_skew( rbtdef.l[i]) ).row_join(sympy.eye(3)*rbtdef.m[i]) )

  M = sympy.zeros((rbtdef.dof,rbtdef.dof))

  rbtdeftmp = copy.deepcopy( rbtdef )
  rbtdeftmp.gravity = sympy.zeros((3,1))
  rbtdeftmp.dq = sympy.zeros((rbtdeftmp.dof,1))

  for i in range( M.rows ):
    rbtdeftmp.ddq = sympy.zeros((rbtdeftmp.dof,1))
    rbtdeftmp.ddq[i] = 1
    geomtmp = geometry.Geometry(rbtdeftmp)

    V, dV = _forward_rne( rbtdeftmp, geomtmp, ifunc )
    Mcoli = _backward_rne( rbtdeftmp, geomtmp, Llm, V, dV, ifunc )

    # It's done like this since M is symmetric:
    M[:,i] = ( M[i,:i].T ) .col_join( Mcoli[i:,:] )

  return M



def _gen_fricterm(rbtdef, ifunc=None):
  '''Generate friction forces (simple Coulomb and viscouse model).'''
  if not ifunc: ifunc = _id
  fric = sympy.zeros((rbtdef.dof,1))
  for i in range( rbtdef.dof ):
    fric[i] = ifunc( rbtdef.fv[i] * rbtdef.dq[i] + rbtdef.fc[i] * sympy.sign(rbtdef.dq[i]) )
  return fric



def _find_dyn_parm_deps( dof, parm_num, regressor_func ):
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
  
  def __init__(self, rbtdef, geom):
    
    self.rbtdef = rbtdef
    self.geom = geom
    self.dof = rbtdef.dof
    
    self.dynparms = sympy.Matrix( rbtdef.dynparms() )
    self.n_dynparms =  len( self.dynparms )
    
    self.delta = self.dynparms
    self.n_delta = self.n_dynparms
    
  def gen_tau(self, ifunc=None):
    self.tau = _gen_tau_rne(self.rbtdef, self.geom, ifunc)
    
  def gen_gravterm(self, ifunc=None):
    self.g = _gen_gravterm_rne(self.rbtdef, self.geom, ifunc)
    
  def gen_ccfterm(self, ifunc=None):
    self.c = _gen_ccfterm_rne(self.rbtdef, self.geom, ifunc)
    
  def gen_massmatrix(self, ifunc=None):
    self.M = _gen_massmatrix_rne(self.rbtdef, self.geom, ifunc)

  def gen_regressor(self, ifunc=None):
    self.regressor = _gen_regressor_rne(self.rbtdef, self.geom, ifunc)
    self.H = self.regressor
    
  def gen_fricterm(self, ifunc=None):
    if self.rbtdef.frictionmodel == 'simple':
      self.f = _gen_fricterm(self.rbtdef, ifunc)
    else:
      self.f = sympy.zeros(self.dof,1)

  def calc_base_parms(self, regressor_func=None):

    if regressor_func == None:
      se = symcode.subexprs.Subexprs(mode='unique_ops')
      regressor = _gen_regressor_rne(self.rbtdef, self.geom, ifunc=se.collect)
      code = (se.subexprs, sympy.flatten(regressor))
      func_def_regressor = symcode.generation.code_to_func('python', code, 'local_regressor_func', ['q','dq','ddq'], [('q'+str(i+1), 'q['+str(i)+']') for i in range(self.dof)])
      global sin, cos, sign
      sin = numpy.sin
      cos = numpy.cos
      sign = numpy.sign
      exec(func_def_regressor)
      regressor_func = local_regressor_func
    
    Pb, Pd, Kd = _find_dyn_parm_deps( self.dof, self.n_delta, regressor_func )

    self.Pb = sympy.Matrix(Pb).applyfunc(lambda x: x.nsimplify())
    self.Pd = sympy.Matrix(Pd).applyfunc(lambda x: x.nsimplify())
    self.Kd = sympy.Matrix(Kd).applyfunc(lambda x: x.nsimplify())

    self.base_idxs = ( numpy.matrix([[i for i in range(self.n_delta)]]) * numpy.matrix(Pb) ).astype(float).astype(int).tolist()[0]

    self.baseparms = ( self.Pb.T + self.Kd * self.Pd.T ) * self.delta
    self.n_base = len( self.baseparms )
    self.base_regressor = self.regressor * self.Pb
    
    self.beta = self.baseparms
    self.n_beta = self.n_base
    self.Hb = self.base_regressor
    
  def gen_all(self, ifunc=None):
    self.gen_tau(ifunc)
    self.gen_gravterm(ifunc)
    self.gen_ccfterm(ifunc)
    self.gen_massmatrix(ifunc)
    self.gen_regressor(ifunc)
    self.calc_base_parms()