### Put this script in the same folder as the sympybotics folder OR set
SYMPYBOTICS_PATH = '/path/to/package/'
# to the path where the SymPyBotics package is (the package is the main folder, inside where __init__.py is)
### To run it, do 'python example.py'

import sys
sys.path += [SYMPYBOTICS_PATH]

import sympybotics

### numpy and sympy libraries are needed
import numpy
import sympy



pi = sympy.pi
q = sympybotics.Robot.q

### First, create a robot object (example of a 2 DOF robot):
rbt = sympybotics.Robot('Example Robot', ### robot name
                        [ (-pi/2 ,0 , 0, q ), ( pi/2 , 0, 0, q )],  ### list of tuples with Denavit-Hartenberg parameters (alpha, a, d, theta) in standard notation
                        'examprobot' ### robot short name (optional)
                       )



### Inside Robot object there are several variables (numerical and symbolic) which represent some quantities:
#  
#  rbt.dof is the DOF
#  rbt.q is the symbolic joint vector
#  rbt.dq and rbt.ddq are first and second derivatives of q
#  rbt.gravity is the gravity acceleration vector (which defaults to [0.0,0.0,-9.81] and
#                                                 can be set with rbt.gravity = sympy.Matrix([...,...,...])
#
#  dynamic symbols:
#  rbt.I is a list of symbolic link inertia tensors
#  rbt.r is a list of symbolic link center of mass vectors
#  rbt.m is a list of symbolis of links masses
#  rbt.L is a list of symbolic inertia tensors about link frames
#  rbt.l is a list of symbolic first moments of inertia (m*r)
#  
#  rbt.dynparms() is a function which return the symbolic dynamic parameters vector
#                 - be aware that L_*xy, L_*xz and L_*yz are the off-diagonal elements
#                   of tensor matrices L_* rather than the "products of inertia" (which
#                   are symmetric)



### Generate geometric model (kinematic model) equations:
rbt.gen_geometric_model()



### Geometric model will be within the rbt.geom object
#
#  rbt.geom.Tdhi is a list with symbolic transformations between consecutive links
#  rbt.geom.Ti is a list with symbolic transformations from links to base
#  rbt.geom.Rdhi and rbt.geom.Ri are lists with respective rotation matrices
#  rbt.geom.pdhi and rbt.geom.pi are lists with respective translation vectors
#
#  For instance, rbt.geom.Ti[2] gives end-effector to base transformation:
#
#  [cos(q1)*cos(q2), -sin(q1), sin(q2)*cos(q1), 0]
#  [sin(q1)*cos(q2),  cos(q1), sin(q1)*sin(q2), 0]
#  [       -sin(q2),        0,         cos(q2), 0]
#  [              0,        0,               0, 1]



### Generate kinemtic model (differential kinematic model) equations:
rbt.gen_kinematic_model()



### Kinematic model will be within the rbt.kinem object
#
#  rbt.kinem.Ji is a list with symbolic Jacobians
#
#  For instance, rbt.geom.Ji[2] gives the end-effector Jacobian:
#
#  [0,        0]
#  [0,        0]
#  [0,        0]
#  [0, -sin(q1)]
#  [0,  cos(q1)]
#  [1,        0]



### Generate dynamic model terms (for higher DOF robots this step can take several minutes):
### This step generates symbolic code rather than symbolic equations since equations became too big
#    to be symbolically treated, although it is also possible to generate such equations.
rbt.dyn = sympybotics.dynamic.Dyn(rbt, usefricdyn=False) # (in this example, such model is saved as rbt.dyn variable)



### Dynamic model code will be within the rbt.dyn object
#
#  rbt.dyn.M_code has the code to compute the mass matrix M(q)
#  rbt.dyn.c_code has the code to compute the Coriolis and centripetal forces term c(q,dq)
#  rbt.dyn.g_code has the code to compute the gravity term g(q)
#  rbt.dyn.tau_code has the code to compute the generalised force (tau) given by M(q)*ddq + c(q,dq) + g(q)
#  rbt.dyn.regressor_code has the code to compute the regressor matrix,
#             i.e., H(q,dq,ddq) in H*delta = tau, where delta is the dynamic parameter vector (as given by rbt.dynparms())
#
### A code object is a pair of two elements:
#   - the first element is a list of pairs, each containing a symbolic intermediate variable and its corresponding symbolic expression
#   - the second element is a list of output expressions (which forms the elements of a vector or a matrix) dependent on the intermediate variables



### To generate C functions from the code do (it's also possible to generate python function):
tau_str = sympybotics.codegen_robot.dyn_code_to_func( 'c', rbt.dyn.tau_code, 'tau', 2, rbt.dof, rbt.dynparms() )
regressor_str = sympybotics.codegen_robot.dyn_code_to_func( 'c', rbt.dyn.regressor_code, 'regressor', 2, rbt.dof )
M_str = sympybotics.codegen_robot.dyn_code_to_func( 'c', rbt.dyn.M_code, 'M', 0, rbt.dof, rbt.dynparms() )
g_str = sympybotics.codegen_robot.dyn_code_to_func( 'c', rbt.dyn.c_code, 'c', 1, rbt.dof, rbt.dynparms() )
c_str = sympybotics.codegen_robot.dyn_code_to_func( 'c', rbt.dyn.g_code, 'g', 0, rbt.dof, rbt.dynparms() )



### For example, M_str will contain a C function to compute the mass matrix
#   where parameter parms must be an array with dynamic parameters values (ordered according to rbt.dynparms())
#   parameter q must be an array with joint positions,
#   and parameter M_out must be an array where the mass matrix will be written in row major order
print('\n\n\n')
print(tau_str)
print('\n\n\n')
print(regressor_str)
print('\n\n\n')
print(M_str)
print('\n\n\n')
print(g_str)
print('\n\n\n')
print(c_str)




