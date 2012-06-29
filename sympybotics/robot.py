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


class Robot(object):
  
  def __init__(self,dof,name,shortname=None):
    
    self.dof = int(dof)
    self.name = str(name)
    if shortname != None :
      self.shortname = str(shortname).replace(' ','_').replace('.','_')
    else :
      self.shortname = name.replace(' ','_').replace('.','_')
    
    
    alpha , a , d , theta = sympy.symbols('alpha,a,d,theta',real=True)
    default_dh_symbols = [ alpha , a , d , theta  ]

    cos = sympy.cos
    sin = sympy.sin
    default_dh_transfmat = sympy.Matrix( [[ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta) ],
                                         [ sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta) ],
                                         [          0,             sin(alpha),             cos(alpha),            d ],
                                         [          0,                      0,                      0,            1 ]] )
    
    #g_a = sympy.symbols('g_a',real=True)
    g_a = 9.81
    self.gravity = sympy.Matrix([[0.0],[0.0],[-g_a]])
    
    self.dh_symbols = default_dh_symbols
    self.dh_transfmat = default_dh_transfmat
    
    self.dyn_parms_order = 'Khalil'
    
    self._gen_symbols()
  
  
  def __repr__(self) :
    return 'Robot instance: ' + self.name
  
  
  def _gen_symbols( self ) :
    
    subs_dict = collections.OrderedDict
    
    def new_sym( name ):
      return sympy.Symbol(name, real=True)

    def sym_skew(v):
      return sympy.Matrix( [ [     0, -v[2],  v[1] ],
                             [  v[2],     0, -v[0] ],
                             [ -v[1],  v[0],     0 ] ] )

    dof = self.dof
    
    #q = self.q = [ new_sym('q_'+str(i+1)) for i in range(self.dof) ]
    #dq = self.dq = [ new_sym('dq_'+str(i+1)) for i in range(self.dof) ]
    #ddq = self.ddq = [ new_sym('ddq_'+str(i+1)) for i in range(self.dof) ]
    q = self.q = sympy.Matrix( [ [new_sym('q'+str(i+1))] for i in range(self.dof) ] )
    dq = self.dq = sympy.Matrix( [ [new_sym('dq'+str(i+1))] for i in range(self.dof) ] )
    ddq = self.ddq = sympy.Matrix( [ [new_sym('ddq'+str(i+1))] for i in range(self.dof) ] )
    
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
      
      m[i] = new_sym('m_'+str(i+1))
      
      r[i] = sympy.Matrix( [ new_sym( 'r_'+str(i+1)+dim ) for dim in ['x','y','z'] ] )
      
      mr[i] = m[i] * r[i]
      l[i] = sympy.Matrix( [ new_sym( 'l_'+str(i+1)+dim ) for dim in ['x','y','z'] ] )
      
      Is[i] = [ new_sym( 'I_'+str(i+1)+elem ) for elem in ['xx','xy','xz', 'yy', 'yz', 'zz'] ]
      Ls[i] = [ new_sym( 'L_'+str(i+1)+elem ) for elem in ['xx','xy','xz', 'yy', 'yz', 'zz'] ]

      I[i] = sympy.Matrix( [ [ Is[i][0], Is[i][1], Is[i][2] ],
                             [ Is[i][1], Is[i][3], Is[i][4] ],
                             [ Is[i][2], Is[i][4], Is[i][5] ] ] )

      L[i] = sympy.Matrix( [ [ Ls[i][0], Ls[i][1], Ls[i][2] ],
                             [ Ls[i][1], Ls[i][3], Ls[i][4] ],
                             [ Ls[i][2], Ls[i][4], Ls[i][5] ] ] )
      
      fv[i] = new_sym( 'fv_'+str(i+1))
      fc[i] = new_sym( 'fc_'+str(i+1))

      
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
    
    
    
  def set_dh_parms( self, dh_parms_list ):
    
    if len( dh_parms_list ) != self.dof:
      raise Exception('Robot.set_geometry(): provided number of links differ from robot dof (%d vs %d).' % ( len( dh_parms_list ), self.dof) )
      
    self.dh_parms = []
    
    for i in range( self.dof ):
      if len( dh_parms_list[i] ) != 4:
        raise Exception('Robot.set_geometry(): wrong number of Denavit-Hartenberg parameters (must be 4 per link).' )
      self.dh_parms.append( dh_parms_list[i] )
      
    return self
      

  def dynparms( self, parm_order = None, usefricdyn=False ):

    if not parm_order: parm_order = self.dyn_parms_order
    parm_order = parm_order.lower()
 
    parms = []
    for i in range( 0, self.dof ):

      if parm_order == 'khalil' or parm_order == 'tensor first': # Lxx Lxy Lxz Lyy Lyz Lzz lx ly lz m
        parms += self.Ls[i]
        parms += self.l[i].mat
        parms += [ self.m[i] ]

      elif parm_order == 'siciliano' or parm_order == 'mass first': # m lx ly lz Lxx Lxy Lxz Lyy Lyz Lzz
        parms += [ self.m[i] ]
        parms += self.l[i].mat
        parms += self.Ls[i]

      else:
          raise Exception('Robot.Parms(): dynamic parameters order \'' +parm_order+ '\' not know.')

      if usefricdyn: parms += [ self.fv[i], self.fc[i] ]

    return parms



  
  
  
  def gen_geometric_model(self):

    from . import geomkinem

    self.geom = geomkinem.gen_geometric_model(self)
    
    self.Tdhi = self.geom.Tdhi
    self.Tdhi_inv = self.geom.Tdhi_inv
    self.Rdhi = self.geom.Rdhi
    self.pdhi = self.geom.pdhi
    self.Ti = self.geom.Ti
    self.Ri = self.geom.Ri
    self.pi = self.geom.pi
    self.zi = self.geom.zi
    self.Si = self.geom.Si
    
    return self


  def Ji( self, link=None ):
    if link == None : link = self.dof
    return self.Jpi[link].col_join( self.Joi[link] )
  

  def gen_kinematic_model( self ):

    from . import geomkinem
    
    self.kinem = geomkinem.gen_kinematic_model(self,self.geom)
    
    self.Jpi = self.kinem.Jpi
    self.Joi = self.kinem.Joi
    self.Jcpi = self.kinem.Jcpi
    self.Jcoi = self.kinem.Jcoi

    return self



