

import sympy



def _new_sym( name ):
  return sympy.symbols(name, real=True)

def _elements_to_tensor(elems):
  return sympy.Matrix( [ [ elems[0], elems[1], elems[2] ],
                         [ elems[1], elems[3], elems[4] ],
                         [ elems[2], elems[4], elems[5] ] ] )

class _elementslist_to_tensorlist():
  def __init__(self,elementslist):
    self.l = elementslist
  def __getitem__(self,i):
    return sympy.Matrix( [ [ self.l[i][0], self.l[i][1], self.l[i][2] ],
                           [ self.l[i][1], self.l[i][3], self.l[i][4] ],
                           [ self.l[i][2], self.l[i][4], self.l[i][5] ] ] )

def _sym_skew(v):
    return sympy.Matrix( [ [     0, -v[2],  v[1] ],
                           [  v[2],     0, -v[0] ],
                           [ -v[1],  v[0],     0 ] ] )

_joint_symb = _new_sym('q')

def _joint_i_symb(i):
  return _new_sym( 'q'+str(i) )

class RobotDef(object):
  """Class that generates and holds robot definitions and symbols."""

  q = _joint_symb

  _dh_alpha , _dh_a , _dh_d , _dh_theta = _new_sym('alpha,a,d,theta')
  _default_dh_symbols = ( _dh_alpha , _dh_a , _dh_d , _dh_theta  )

  _cos = sympy.cos
  _sin = sympy.sin
  _default_dh_transfmat = sympy.Matrix( [ 
    [ _cos(_dh_theta), -_sin(_dh_theta)*_cos(_dh_alpha),  _sin(_dh_theta)*_sin(_dh_alpha), _dh_a*_cos(_dh_theta) ],
    [ _sin(_dh_theta),  _cos(_dh_theta)*_cos(_dh_alpha), -_cos(_dh_theta)*_sin(_dh_alpha), _dh_a*_sin(_dh_theta) ],
    [          0,             _sin(_dh_alpha),             _cos(_dh_alpha),            _dh_d ],
    [          0,                      0,                      0,            1 ]
    ] )
    
  
  def __init__(self,name,dh_parms,shortname=None):
    """
    Create RobotDef instance with data structures for robot geometry and symbols of robot dynamics.
    """
    
    self.dof = len(dh_parms)
    
    self.name = str(name)
    if shortname != None :
      self.shortname = str(shortname).replace(' ','_').replace('.','_')
    else :
      self.shortname = ''.join(c for c in name.replace(' ','_').replace('.','_') if c.isalnum() or c == '_')
    
    self.dh_symbols = RobotDef._default_dh_symbols
    self.dh_transfmat = RobotDef._default_dh_transfmat
    
    self.dyn_parms_order = 'Khalil'
    
    self.frictionmodel = None # can be None or 'simple'

    #g_a = sympy.symbols('g_a',real=True)
    g_a = 9.81
    self.gravity = sympy.Matrix([[0.0],[0.0],[-g_a]])
    
    
    self._gen_symbols()

    self._set_dh_parms(dh_parms)
  
  def __setattr__(self, name, value):
    object.__setattr__(self, name, value)
    if name == 'Le':
      object.__setattr__(self, 'L', _elementslist_to_tensorlist(value))
    elif name == 'Ie':
      object.__setattr__(self, 'I', _elementslist_to_tensorlist(value))
  
  def __repr__(self) :
    return 'RobotDef instance: ' + self.name

  
  def _gen_symbols( self ) :
    """Generate robot dynamic symbols and populates RobotDef instance with them. (internal function)"""
     

    dof = self.dof
    
    #q = self.q = [ _new_sym('q_'+str(i+1)) for i in range(self.dof) ]
    #dq = self.dq = [ _new_sym('dq_'+str(i+1)) for i in range(self.dof) ]
    #ddq = self.ddq = [ _new_sym('ddq_'+str(i+1)) for i in range(self.dof) ]
    q = self.q = sympy.Matrix( [ [_new_sym('q'+str(i+1))] for i in range(self.dof) ] )
    dq = self.dq = sympy.Matrix( [ [_new_sym('dq'+str(i+1))] for i in range(self.dof) ] )
    ddq = self.ddq = sympy.Matrix( [ [_new_sym('ddq'+str(i+1))] for i in range(self.dof) ] )
    
    m = self.m = list( range( self.dof ) )
    l = self.l = list( range( self.dof ) )
    Le = self.Le = list( range( self.dof ) )

    L = self.L

    r = self.r = list( range( self.dof ) )
    Ie = self.Ie = list( range( self.dof ) )
    
    I = self.I
    
    fv = self.fv = list( range( self.dof ) )
    fc = self.fc = list( range( self.dof ) )
    
    I_funcof_L = self.I_funcof_L = list( range( self.dof ) )
    L_funcof_I = self.L_funcof_I = list( range( self.dof ) )
    
    dict_I2Lexp = self.dict_I2Lexp = dict()
    dict_L2Iexp = self.dict_L2Iexp = dict()
    dict_l2mr = self.dict_l2mr = dict()
    dict_r2lm = self.dict_r2lm = dict()
    
    for i in range( dof ):
      
      m[i] = _new_sym('m_'+str(i+1))
      l[i] = sympy.Matrix( [ _new_sym( 'l_'+str(i+1)+dim ) for dim in ['x','y','z'] ] )
      Le[i] = [ _new_sym( 'L_'+str(i+1)+elem ) for elem in ['xx','xy','xz', 'yy', 'yz', 'zz'] ]
      
      r[i] = sympy.Matrix( [ _new_sym( 'r_'+str(i+1)+dim ) for dim in ['x','y','z'] ] )      
      Ie[i] = [ _new_sym( 'I_'+str(i+1)+elem ) for elem in ['xx','xy','xz', 'yy', 'yz', 'zz'] ]
      
      fv[i] = _new_sym( 'fv_'+str(i+1))
      fc[i] = _new_sym( 'fc_'+str(i+1))

      L_funcof_I[i] = I[i] + m[i] * _sym_skew(r[i]).T * _sym_skew(r[i])
      I_funcof_L[i] = L[i] - m[i] * _sym_skew(r[i]).T * _sym_skew(r[i])

      for elem,exprss in enumerate(I_funcof_L[i]):
        dict_I2Lexp[ I[i][elem] ] = exprss

      for elem,exprss in enumerate(L_funcof_I[i]):
        dict_L2Iexp[ L[i][elem] ] = exprss

      for elem in range(3):
        dict_l2mr[ l[i][elem] ] = m[i] * r[i][elem]
        dict_r2lm[ r[i][elem] ] = l[i][elem] / m[i]

      self.latex_symbols = {}
      for i in range(self.dof):
        self.latex_symbols[ self.q[i] ] = r'q_{'+str(i+1)+'}'
        self.latex_symbols[ self.dq[i] ] = r'\dot{q}_{'+str(i+1)+'}'
        self.latex_symbols[ self.ddq[i] ] = r'\ddot{q}_{'+str(i+1)+'}'

    return self
    
    
    
  def _set_dh_parms( self, dh_parms_list ):
    """
    Define the RobotDef geometry using Denavit-Hartenberg notation.
    """
    
    if len( dh_parms_list ) != self.dof:
      raise Exception('RobotDef.set_geometry(): provided number of links differ from robot dof (%d vs %d).' % ( len( dh_parms_list ), self.dof) )
      
    self.dh_parms = []

    self.links_sigma = [0]*self.dof
    theta_index = self.dh_symbols.index(sympy.Symbol('theta',real=True))
    d_index = self.dh_symbols.index(sympy.Symbol('d',real=True))
    
    for i in range( self.dof ):
      
      if len( dh_parms_list[i] ) != 4:
        raise Exception('RobotDef.set_dh_parms: wrong number of Denavit-Hartenberg parameters (must be 4 per link).' )

      for j,p in enumerate(dh_parms_list[i]):

        for v in sympy.sympify(p).free_symbols:
          
          v=str(v)
          if v[0] == str(_joint_symb):
            
            if len(v) > 1:
              try:
                num = int(v[1:])
              except:
                num = 1
              if num <= 0 or num > self.dof:
                raise Exception("RobotDef.set_dh_parms: Joint position symbol \'%s\' out of robot joint range (from 1 to %d)!" % (v,self.dof))

            else:
              temp = list(dh_parms_list[i])
              temp[j] = sympy.sympify( temp[j] ).subs( {_joint_symb:_joint_i_symb(i+1)} )
              dh_parms_list[i] = tuple(temp)
            
      self.dh_parms.append( dh_parms_list[i] )
      
      try:
        if dh_parms_list[i][ theta_index ].has( self.q[i] ):
          self.links_sigma[i] = 0
          # print 'joint',i+1,'is revolute'
      except: pass
      try:
        if dh_parms_list[i][ d_index ].has( self.q[il] ):
          self.links_sigma[i] = 1
          # print 'joint',il+1,'is prismatic'
      except: pass
      
    return self
      

  def dynparms( self, parm_order = None ):
    """Return list of RobotDef symbolic dynamic parameters."""

    if not parm_order: parm_order = self.dyn_parms_order
    parm_order = parm_order.lower()
 
    parms = []
    for i in range( 0, self.dof ):

      if parm_order == 'khalil' or parm_order == 'tensor first': # Lxx Lxy Lxz Lyy Lyz Lzz lx ly lz m
        parms += self.Le[i]
        parms += sympy.flatten(self.l[i])
        parms += [ self.m[i] ]

      elif parm_order == 'siciliano' or parm_order == 'mass first': # m lx ly lz Lxx Lxy Lxz Lyy Lyz Lzz
        parms += [ self.m[i] ]
        parms += sympy.flatten(self.l[i])
        parms += self.Le[i]

      else:
          raise Exception('RobotDef.Parms(): dynamic parameters order \'' +parm_order+ '\' not know.')

      if self.frictionmodel == 'simple': parms += [ self.fv[i], self.fc[i] ]

    return parms


