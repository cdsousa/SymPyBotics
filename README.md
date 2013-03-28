SymPyBotics
===========

Symbolic Robotics Toolbox using Python and SymPy

Example
-------

Definition of a 2 DOF example robot:

~~~~~~~~~~~Python
>>> import sympy
>>> import sympybotics
>>> pi = sympy.pi
>>> q = sympybotics.RobotDef.q
>>> rbtdef = sympybotics.RobotDef('Example Robot', # robot name
                              [ (-pi/2 ,0 , 0, q+pi/2 ),   # list of tuples with standard Denavit-Hartenberg parameters 
                                ( pi/2 , 0, 0, q-pi/2 ) ], # (alpha, a, d, theta)
                              'examprobot' # robot short name (optional)
                             )
~~~~~~~~~~~

~~~~~~~~~~~Python
>>> print(rbtdef.dynparms())
[L_1xx, L_1xy, L_1xz, L_1yy, L_1yz, L_1zz, l_1x, l_1y, l_1z, m_1, L_2xx, L_2xy, L_2xz, L_2yy, L_2yz, L_2zz, l_2x, l_2y, l_2z, m_2]
~~~~~~~~~~~

Generation of geometric, kinematic and dynamic models:

~~~~~~~~~~~Python
>>> rbt = sympybotics.RobotDynCode(rbtdef, codecollectmode='unique_ops')
generating geometric model
generating kinematic model
generating tau code
generating gravity term code
generating coriolis term code
generating inertia matrix code
generating regressor matrix code
done
~~~~~~~~~~~

~~~~~~~~~~~Python
>>> print(rbt.geom.T[-1])
[-sin(q1)*sin(q2), -cos(q1),  sin(q1)*cos(q2), 0]
[ sin(q2)*cos(q1), -sin(q1), -cos(q1)*cos(q2), 0]
[         cos(q2),        0,          sin(q2), 0]
[               0,        0,                0, 1]
~~~~~~~~~~~

~~~~~~~~~~~Python
>>> rbt.kin.J[-1]
[0,        0]
[0,        0]
[0,        0]
[0, -cos(q1)]
[0, -sin(q1)]
[1,        0]
~~~~~~~~~~~

Dynamic code optimization:

~~~~~~~~~~~Python
>>> rbt.optimize_code(mode='light')
optimizing tau_code
optimizing g_code
optimizing c_code
optimizing M_code
optimizing H_code
done
~~~~~~~~~~~

~~~~~~~~~~~Python
>>> tau_str = sympybotics.robotcodegen.dyn_code_to_func('c', rbt.tau_code, 'tau', 2, rbt.dof, rbtdef.dynparms())
>>> print(tau_str)
void tau( double* tau_out, const double* parms, const double* q, const double* dq, const double* ddq )
{
  double tmp0 = -dq[0];
  double tmp1 = -ddq[0];
  double tmp2 = cos(q[1]);
  double tmp3 = -tmp0*tmp2;
  double tmp4 = sin(q[1]);
  double tmp5 = -tmp0*tmp4;
  double tmp8 = dq[1]*tmp0*tmp4 - tmp1*tmp2;
  double tmp11 = -dq[1]*tmp0*tmp2 - tmp1*tmp4;
  double tmp12 = 9.81*tmp2;
  double tmp13 = 9.81*tmp4;
  double tmp17 = parms[12]*tmp3 + parms[14]*dq[1] + parms[15]*tmp5;
  double tmp24 = parms[11]*tmp3 + parms[13]*dq[1] + parms[14]*tmp5;
  double tmp38 = parms[10]*tmp3 + parms[11]*dq[1] + parms[12]*tmp5;

  tau_out[0] = -parms[3]*tmp1 + tmp2*(parms[10]*tmp8 + parms[11]*ddq[1] + parms[12]*tmp11 + dq[1]*tmp17 + parms[17]*tmp13 - tmp24*tmp5) + tmp4*(parms[12]*tmp8 + parms[14]*ddq[1] + parms[15]*tmp11 - dq[1]*tmp38 - parms[17]*tmp12 + tmp24*tmp3);
  tau_out[1] = parms[11]*tmp8 + parms[13]*ddq[1] + parms[14]*tmp11 - parms[16]*tmp13 + parms[18]*tmp12 - tmp17*tmp3 + tmp38*tmp5;

  return;
}
~~~~~~~~~~~

Dynamic base parameters:

~~~~~~~~~~~
>>> rbt.calc_base_parms()
calculating base parameters and regressor code
done
>>> print(rbt.dyn.baseparms)
[L_1yy + L_2zz]
[L_2xx - L_2zz]
[        L_2xy]
[        L_2xz]
[        L_2yy]
[        L_2yz]
[         l_2x]
[         l_2z]
~~~~~~~~~~~

