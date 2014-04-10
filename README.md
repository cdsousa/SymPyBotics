SymPyBotics
===========

Symbolic Framework for Modeling and Identification of Robot Dynamics

Uses [Sympy](http://sympy.org) and [Numpy](http://www.numpy.org/) libraries.

[![Build Status](https://travis-ci.org/cdsousa/SymPyBotics.png?branch=master)](https://travis-ci.org/cdsousa/SymPyBotics)

Example
-------

Definition of a 2 DOF example robot:

```Python
>>> import sympy
>>> import sympybotics
>>> rbtdef = sympybotics.RobotDef('Example Robot', # robot name
...                               [('-pi/2', 0, 0, 'q+pi/2'),  # list of tuples with Denavit-Hartenberg parameters
...                                ( 'pi/2', 0, 0, 'q-pi/2')], # (alpha, a, d, theta)
...                               dh_convention='standard' # either 'standard' or 'modified'
...                              )
>>> rbtdef.frictionmodel = {'Coulomb', 'viscous'} # options are None or a combination of 'Coulomb', 'viscous' and 'offset'

```

```Python
>>> rbtdef.dynparms()
[L_1xx, L_1xy, L_1xz, L_1yy, L_1yz, L_1zz, l_1x, l_1y, l_1z, m_1, fv_1, fc_1, L_2xx, L_2xy, L_2xz, L_2yy, L_2yz, L_2zz, l_2x, l_2y, l_2z, m_2, fv_2, fc_2]

```

Generation of geometric, kinematic and dynamic models:

```Python
>>> rbt = sympybotics.RobotDynCode(rbtdef, verbose=True)
generating geometric model
generating kinematic model
generating inverse dynamics code
generating gravity term code
generating coriolis term code
generating inertia matrix code
generating regressor matrix code
generating friction term code
done

```

```Python
>>> rbt.geo.T[-1]
Matrix([
[-sin(q1)*sin(q2), -cos(q1),  sin(q1)*cos(q2), 0],
[ sin(q2)*cos(q1), -sin(q1), -cos(q1)*cos(q2), 0],
[         cos(q2),        0,          sin(q2), 0],
[               0,        0,                0, 1]])

```

```Python
>>> rbt.kin.J[-1]
Matrix([
[0,        0],
[0,        0],
[0,        0],
[0, -cos(q1)],
[0, -sin(q1)],
[1,        0]])

```

C function generation:

```Python
>>> tau_str = sympybotics.robotcodegen.robot_code_to_func('C', rbt.invdyn_code, 'tau_out', 'tau', rbtdef)


```
Doing `print(tau_str)`, function code will be output:

```C
void tau( double* tau_out, const double* parms, const double* q, const double* dq, const double* ddq )
{
  double x0 = sin(q[1]);
  double x1 = -dq[0];
  double x2 = -x1;
  double x3 = x0*x2;
  double x4 = cos(q[1]);
  double x5 = x2*x4;
  double x6 = parms[13]*x5 + parms[15]*dq[1] + parms[16]*x3;
  double x7 = parms[14]*x5 + parms[16]*dq[1] + parms[17]*x3;
  double x8 = -ddq[0];
  double x9 = -x4;
  double x10 = dq[1]*x1;
  double x11 = x0*x10 + x8*x9;
  double x12 = -x0*x8 - x10*x4;
  double x13 = 9.81*x0;
  double x14 = 9.81*x4;
  double x15 = parms[12]*x5 + parms[13]*dq[1] + parms[14]*x3;

  tau_out[0] = -parms[3]*x8 + x0*(parms[14]*x11 + parms[16]*ddq[1] + parms[17]*x12 - dq[1]*x15 - parms[19]*x14 + x5*x6) - x9*(parms[12]*x11 + parms[13]*ddq[1] + parms[14]*x12 + dq[1]*x7 + parms[19]*x13 - x3*x6);
  tau_out[1] = parms[13]*x11 + parms[15]*ddq[1] + parms[16]*x12 - parms[18]*x13 + parms[20]*x14 + x15*x3 - x5*x7;

  return;
}
```

Dynamic base parameters:

```Python
>>> rbt.calc_base_parms()
>>> rbt.dyn.baseparms
Matrix([
[L_1yy + L_2zz],
[         fv_1],
[         fc_1],
[L_2xx - L_2zz],
[        L_2xy],
[        L_2xz],
[        L_2yy],
[        L_2yz],
[         l_2x],
[         l_2z],
[         fv_2],
[         fc_2]])

```

Author
------

[Cristóvão Duarte Sousa](https://github.com/cdsousa)

Install
-------

From git source:

    git clone https://github.com/cdsousa/SymPyBotics.git
    cd sympybotics
    python setup.py install

License
-------

New BSD license. See [License File](LICENSE.txt)
