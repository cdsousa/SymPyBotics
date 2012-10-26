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


import collections
import sympy



def _new_sym( name ):
  return sympy.symbols(name, real=True)


class Robot(object):
  """Class that generates and holds robot geometric and kinematic information, and dynamic parameter symbols."""


  _joint_pos_symb = _new_sym('q')

  q = _joint_pos_symb

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
    Create Robot instance with data structures for robot geometry and symbols of robot dynamics.
    """
    
    self.dof = len(dh_parms)
    
    self.name = str(name)
    if shortname != None :
      self.shortname = str(shortname).replace(' ','_').replace('.','_')
    else :
      self.shortname = ''.join(c for c in name.replace(' ','_').replace('.','_') if c.isalnum() or c == '_')
    
    self.dh_symbols = Robot._default_dh_symbols
    self.dh_transfmat = Robot._default_dh_transfmat
    
    self.dyn_parms_order = 'Khalil'

    #g_a = sympy.symbols('g_a',real=True)
    g_a = 9.81
    self.gravity = sympy.Matrix([[0.0],[0.0],[-g_a]])
    
    
    self._gen_symbols()

    self._set_dh_parms(dh_parms)
  

  
  def __repr__(self) :
    return 'Robot instance: ' + self.name
  
  
  def _gen_symbols( self ) :
    """Generate robot dynamic symbols and populates Robot instance with them. (internal function)"""
    
    subs_dict = collections.OrderedDict
    
    def sym_skew(v):
      return sympy.Matrix( [ [     0, -v[2],  v[1] ],
                             [  v[2],     0, -v[0] ],
                             [ -v[1],  v[0],     0 ] ] )

    dof = self.dof
    
    #q = self.q = [ _new_sym('q_'+str(i+1)) for i in range(self.dof) ]
    #dq = self.dq = [ _new_sym('dq_'+str(i+1)) for i in range(self.dof) ]
    #ddq = self.ddq = [ _new_sym('ddq_'+str(i+1)) for i in range(self.dof) ]
    q = self.q = sympy.Matrix( [ [_new_sym('q'+str(i+1))] for i in range(self.dof) ] )
    dq = self.dq = sympy.Matrix( [ [_new_sym('dq'+str(i+1))] for i in range(self.dof) ] )
    ddq = self.ddq = sympy.Matrix( [ [_new_sym('ddq'+str(i+1))] for i in range(self.dof) ] )
    
    m = self.m = list( range( self.dof ) )

    r = self.r = list( range( self.dof ) )
    
    mr = self.ml = list( range( self.dof ) )
    l = self.l = list( range( self.dof ) )
    
    Is = self.Is = list( range( self.dof ) )
    Ls = self.Ls = list( range( self.dof ) )
    
    I = self.I = list( range( self.dof ) )
    L = self.L = list( range( self.dof ) )
    
    fv = self.fv = list( range( self.dof ) )
    fc = self.fc = list( range( self.dof ) )
    
    I_funcof_L = self.I_funcof_L = list( range( self.dof ) )
    L_funcof_I = self.L_funcof_I = list( range( self.dof ) )
    
    dict_I2Lexp = self.dict_I2Lexp = subs_dict()
    dict_L2Iexp = self.dict_L2Iexp = subs_dict()
    dict_l2mr = self.dict_l2mr = subs_dict()
    dict_r2lm = self.dict_r2lm = subs_dict()


    
    for i in range( dof ):
      
      m[i] = _new_sym('m_'+str(i+1))
      
      r[i] = sympy.Matrix( [ _new_sym( 'r_'+str(i+1)+dim ) for dim in ['x','y','z'] ] )
      
      mr[i] = m[i] * r[i]
      l[i] = sympy.Matrix( [ _new_sym( 'l_'+str(i+1)+dim ) for dim in ['x','y','z'] ] )
      
      Is[i] = [ _new_sym( 'I_'+str(i+1)+elem ) for elem in ['xx','xy','xz', 'yy', 'yz', 'zz'] ]
      Ls[i] = [ _new_sym( 'L_'+str(i+1)+elem ) for elem in ['xx','xy','xz', 'yy', 'yz', 'zz'] ]

      I[i] = sympy.Matrix( [ [ Is[i][0], Is[i][1], Is[i][2] ],
                             [ Is[i][1], Is[i][3], Is[i][4] ],
                             [ Is[i][2], Is[i][4], Is[i][5] ] ] )

      L[i] = sympy.Matrix( [ [ Ls[i][0], Ls[i][1], Ls[i][2] ],
                             [ Ls[i][1], Ls[i][3], Ls[i][4] ],
                             [ Ls[i][2], Ls[i][4], Ls[i][5] ] ] )
      
      fv[i] = _new_sym( 'fv_'+str(i+1))
      fc[i] = _new_sym( 'fc_'+str(i+1))

      
      L_funcof_I[i] = I[i] + m[i] * sym_skew(r[i]).T * sym_skew(r[i])
      I_funcof_L[i] = L[i] - m[i] * sym_skew(r[i]).T * sym_skew(r[i])

      for elem,exprss in enumerate(I_funcof_L[i]):
        dict_I2Lexp[ I[i][elem] ] = exprss

      for elem,exprss in enumerate(L_funcof_I[i]):
        dict_L2Iexp[ L[i][elem] ] = exprss

      for elem in range(3):
        dict_l2mr[ l[i][elem] ] = mr[i][elem]
        dict_r2lm[ r[i][elem] ] = l[i][elem] / m[i]

      self.latex_symbols = {}
      for i in range(self.dof):
        self.latex_symbols[ self.dq[i] ] = r'\dot{q}_'+str(i+1)
        self.latex_symbols[ self.ddq[i] ] = r'\ddot{q}_'+str(i+1)

    return self
    
    
    
  def _set_dh_parms( self, dh_parms_list ):
    """
    Define the Robot geometry using Denavit-Hartenberg notation.
    """
    
    if len( dh_parms_list ) != self.dof:
      raise Exception('Robot.set_geometry(): provided number of links differ from robot dof (%d vs %d).' % ( len( dh_parms_list ), self.dof) )
      
    self.dh_parms = []
    
    for i in range( self.dof ):
      
      if len( dh_parms_list[i] ) != 4:
        raise Exception('Robot.set_dh_parms: wrong number of Denavit-Hartenberg parameters (must be 4 per link).' )

      for j,p in enumerate(dh_parms_list[i]):

        for v in sympy.sympify(p).free_symbols:
          
          v=str(v)
          if v[0] == str(Robot._joint_pos_symb):
            
            if len(v) > 1:
              try:
                num = int(v[1:])
              except:
                num = 1
              if num <= 0 or num > self.dof:
                raise Exception("Robot.set_dh_parms: Joint position symbol \'%s\' out of robot joint range (from 1 to %d)!" % (v,self.dof))

            else:
              temp = list(dh_parms_list[i])
              temp[j] = sympy.sympify( temp[j] ).subs( {Robot._joint_pos_symb:self.q[i]} )
              dh_parms_list[i] = tuple(temp)
            
      self.dh_parms.append( dh_parms_list[i] )
      
    return self
      

  def dynparms( self, parm_order = None, usefricdyn=False ):
    """Return list of Robot symbolic dynamic parameters."""

    if not parm_order: parm_order = self.dyn_parms_order
    parm_order = parm_order.lower()
 
    parms = []
    for i in range( 0, self.dof ):

      if parm_order == 'khalil' or parm_order == 'tensor first': # Lxx Lxy Lxz Lyy Lyz Lzz lx ly lz m
        parms += self.Ls[i]
        parms += sympy.flatten(self.l[i])
        parms += [ self.m[i] ]

      elif parm_order == 'siciliano' or parm_order == 'mass first': # m lx ly lz Lxx Lxy Lxz Lyy Lyz Lzz
        parms += [ self.m[i] ]
        parms += sympy.flatten(self.l[i])
        parms += self.Ls[i]

      else:
          raise Exception('Robot.Parms(): dynamic parameters order \'' +parm_order+ '\' not know.')

      if usefricdyn: parms += [ self.fv[i], self.fc[i] ]

    return parms


  
  def gen_geometric_model(self):
    """Generate symbolic geometric model member object."""
    from . import geomkinem
    self.geom = geomkinem.Geom(self)
    return self


  def gen_kinematic_model( self ):
    """Generate symbolic kinematic model member object."""
    from . import geomkinem
    self.kinem = geomkinem.Kinem(self,self.geom)
    return self



