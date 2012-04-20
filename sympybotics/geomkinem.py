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


import sympy

from . import intermediate


def gen_geometric_model( robot , gen_intervars = False ):
  
  print("Attention to indexes -- dof+1")

  class Geom(object):
    pass

  def inverse_T(T):
    return T[0:3,0:3].transpose().row_join( - T[0:3,0:3].transpose() * T[0:3,3] ).col_join( sympy.zeros((1,3)).row_join(sympy.eye(1)) )
  
  geom = Geom()
  
  ivars = []
  if gen_intervars and not isinstance( gen_intervars , str ): gen_intervars = 'ivargeom_'
  m_intervar_func = intermediate.genfunc_m_intervar( gen_intervars, ivars )
  
  (alpha , a , d , theta) = sympy.symbols('alpha,a,d,theta',real=True)
  
  geom.Tdhi = list(range(robot.dof+1))
  geom.Tdhi[0] = sympy.eye(4 )
  geom.Tdhi_inv = list(range(robot.dof+1))
  geom.Tdhi_inv[0] = sympy.eye(4 )
  geom.Rdhi = list(range(robot.dof+1))
  geom.Rdhi[0] = sympy.eye(3 )
  geom.pdhi = list(range(robot.dof+1))
  geom.pdhi[0] = sympy.zeros((3,1 ))
  
  dh_transfmat_inv = inverse_T(robot.dh_transfmat)
  
  robot.links_sigma = [0]*robot.dof
  
  for l in range(robot.dof):
    
    subs_dict = dict( zip( robot.dh_symbols, robot.dh_parms[l] ) )
    
    geom.Tdhi[l+1] = robot.dh_transfmat.subs(subs_dict)
    geom.Tdhi_inv[l+1] = dh_transfmat_inv.subs(subs_dict)
    geom.Rdhi[l+1] = geom.Tdhi[l+1][0 :3 ,0 :3 ]
    geom.pdhi[l+1] = geom.Tdhi[l+1][0 :3 ,3 ]
    
    try:
      if subs_dict[ theta ].has( robot.q[l,0] ):
        robot.links_sigma[l] = 0
        # print 'joint',l+1,'is revolute'
    except: pass
    try:
      if subs_dict[ d ].has( robot.q[l,0] ):
        robot.links_sigma[l] = 1
        # print 'joint',l+1,'is prismatic'
    except: pass
  
  geom.Ti = list(range(robot.dof+1 ))
  geom.Ti[0] = sympy.eye(4)
  
  for l in range(1 , robot.dof+1 ):
    
    geom.Ti[l] = m_intervar_func( geom.Ti[l-1 ] *  geom.Tdhi[l] )
        
    #Ti[l] = Ti[l].simplify_rational()
    #Ti[l] = trig_reduce(Ti[l])
    #Ti[l] = Ti[l].simplify()
  
  geom.Ri = list(range(0 ,robot.dof+1))
  geom.pi = list(range(0 ,robot.dof+1))
  geom.zi = list(range(0 ,robot.dof+1))
  
  for l in range(0 ,robot.dof+1 ):
    geom.Ri[l] = geom.Ti[l][0 :3 ,0 :3 ]
    geom.pi[l] = geom.Ti[l][0 :3 ,3 ]
    geom.zi[l] = geom.Ri[l][0 :3 ,2 ]
  
  # for screw theory:

  cos = sympy.cos
  sin = sympy.sin
  exp = sympy.exp
  
  Mr = sympy.Matrix([[  1,           0,            0,  a  ],
               [  0,  cos(alpha),  -sin(alpha),  0  ],
               [  0,  sin(alpha),   cos(alpha),  d  ],
               [  0,           0,            0,  1  ]])
  Pr = sympy.Matrix([[  0, -1,  0,  0  ],
               [  1,  0,  0,  0  ],
               [  0,  0,  0,  0  ],
               [  0,  0,  0,  0  ]])
  ### from Frank Park paper:
  #Mp = sympy.Matrix([[  cos(theta),  -sin(theta)*cos(alpha),            0,  a*cos(theta)  ],
               #[  sin(theta),   cos(theta)*cos(alpha),  -sin(alpha),  a*sin(theta)  ],
               #[           0,              sin(alpha),   cos(alpha),             0  ],
               #[           0,                       0,            0,             1  ]])
  ### my own:
  Mp = sympy.Matrix([[  cos(theta),  -sin(theta)*cos(alpha),   sin(theta)*sin(alpha),  a*cos(theta)  ],
               [  sin(theta),   cos(theta)*cos(alpha),  -cos(theta)*sin(alpha),  a*sin(theta)  ],
               [           0,              sin(alpha),              cos(alpha),             0  ],
               [           0,                       0,                       0,             1  ]])
  Pp = sympy.Matrix([[  0,  0,  0,  0  ],
               [  0,  0,  0,  0  ],
               [  0,  0,  0,  1  ],
               [  0,  0,  0,  0  ]])
      
  #if 0:
      #D_exp2trig = { exp(SR.wild(0)) : exp(SR.wild(0).real_part()) * ( cos(SR.wild(0).imag_part()) + sympy.I*sin(SR.wild(0).imag_part()) ) }
      #ePtM_r = exp( Pr*theta ) * Mr
      #ePtM_r = ePtM_r.expand().subs(D_exp2trig).simplify_rational()
      #ePtM_p = exp( Pp*d ) * Mp
      #ePtM_p = ePtM_p.expand().subs(D_exp2trig).simplify_rational()
      #if bool( ePtM_r != robot.dh_transfmat or ePtM_p != robot.dh_transfmat ):
          #raise Exception('gen_geometric_model: interlink transformation does not follows implemented DH formulation')
  
  Sr = ( inverse_T(Mr) * Pr * Mr ).applyfunc(lambda x: sympy.trigsimp(x))
  Sp = ( inverse_T(Mp) * Pp * Mp ).applyfunc(lambda x: sympy.trigsimp(x))
      
      
      
  def sym_se3_unskew(g):
    w = sympy.Matrix( [ g[2,1], g[0,2], g[1,0] ] )
    v = g[0:3,3]
    return w.col_join(v)
        
  geom.Si = list(range(robot.dof+1))
  for l in range(robot.dof):
      if robot.links_sigma[l]: S = Sp
      else: S = Sr
      geom.Si[l+1] = m_intervar_func( sym_se3_unskew( S.subs( dict( zip(robot.dh_symbols, robot.dh_parms[l]) ) ) ) )

  if gen_intervars:
    return ivars,geom
  else:
    return geom





def gen_kinematic_model( robot , geom, gen_intervars = False ):
  
  print("Attention to indexes -- dof+1")
  
  class Kinem(object):
    pass
  
  def sym_skew(v):
    return sympy.Matrix( [ [     0, -v[2],  v[1] ],
    [  v[2],     0, -v[0] ],
    [ -v[1],  v[0],     0 ] ] )
  
  kinem = Kinem()

  ivars = []
  if gen_intervars and not isinstance( gen_intervars , str ): gen_intervars = 'ivarkinm_'
  m_intervar_func = intermediate.genfunc_m_intervar( gen_intervars, ivars )
  
  kinem.Jpi = list(range(0 ,robot.dof+1 ))
  kinem.Jpi[0] = sympy.zeros((3,robot.dof))
  for l in range(1 ,robot.dof+1 ):
    kinem.Jpi[l] = sympy.zeros((3 ,robot.dof))
    for j in range(1 ,l+1 ):
      if robot.links_sigma[j-1]:
        kinem.Jpi[l][0 :3 ,j-1 ] = m_intervar_func( geom.zi[j-1] )
      else:
        kinem.Jpi[l][0 :3 ,j-1 ] = m_intervar_func( geom.zi[j-1].cross( ( geom.pi[l] - geom.pi[j-1] ) ).transpose() )

    #Jpi[i] = Jpi[i].simplify_rational()
    #Jpi[i] = trig_reduce(Jpi[i])
    #Jpi[i] = Jpi[i].simplify()

  kinem.Joi = list(range(0 ,robot.dof+1 ))
  kinem.Joi[0] = sympy.zeros((3,robot.dof))
  for l in range(1 ,robot.dof+1 ):
    kinem.Joi[l] = sympy.zeros((3,robot.dof))
    for j in range(1 ,l+1 ):
      if robot.links_sigma[j-1]:
        kinem.Joi[l][0 :3 ,j-1 ] = sympy.zeros(( 3,1))
      else:
        kinem.Joi[l][0 :3 ,j-1 ] = m_intervar_func( geom.zi[j-1 ] )

    #Joi[i] = Joi[i].simplify_rational()
    #Joi[i] = Joi[i].simplify()

  kinem.Jcpi = list(range(0 ,robot.dof+1 ))
  kinem.Jcoi = kinem.Joi
  for l in range(1 ,robot.dof+1 ):
    kinem.Jcpi[l] = m_intervar_func( kinem.Jpi[l] - sym_skew( geom.Ri[l]*robot.l[l-1] ) * kinem.Joi[l] )
    
  #kinem.Ji = range(0 ,robot.dof+1 )
  #kinem.Ji[0] = sympy.zeros((6,robot.dof))
  #for l in range(1 ,robot.dof+1 ):
     #kinem.Ji[l] = kinem.Jpi[l].stack(kinem.Joi[l])
  
  
  def Ji( robot, link=None ):
    if link == None : link = robot.dof
    return robot.Jpi[link].col_join( robot.Joi[link] )
  
  kinem.Ji = Ji

  if gen_intervars:
    return ivars,kinem
  else:
    return kinem

