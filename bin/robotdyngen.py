#!/usr/bin/env python

from __future__ import print_function

import argparse
import yaml
import sys
import textwrap

from sympy import Matrix
from numpy import sin, cos, sign
from sympybotics import RobotDef, Geometry, Kinematics, Dynamics, robot_code_to_func
from sympybotics.symcode import Subexprs
from sympybotics._compatibility_ import exec_


def fprint(x):
    print(x, end='')
    sys.stdout.flush()


def gen_robot(defs_dict):

    robotdef_d = defs_dict['robot']

    robotdef_kwargs = {}
    robotdef_kwargs['name'] = robotdef_d['name']
    robotdef_kwargs['dh_parms'] = robotdef_d['DH parameters']
    robotdef_kwargs['dh_convention'] = robotdef_d['DH convention']
    robotdef_kwargs['shortname'] = robotdef_d.get('short name', None)

    rbtdef = RobotDef(**robotdef_kwargs)

    if 'extra' in robotdef_d:
        if 'friction model' in robotdef_d['extra']:
            rbtdef.frictionmodel = robotdef_d['extra']['friction model']
        if 'gravity acceleration' in robotdef_d['extra']:
            rbtdef.gravityacc = Matrix(robotdef_d['extra']['gravity acceleration'])

    print(rbtdef.description())

    source = defs_dict['generate']['source']
    outputs = source['outputs']

    outnames_dict = {'transf':'T',
                     'jacob':'J',
                     'invdyn':'tau',
                     'regressor':'H',
                     'inertiamatrix':'M',
                     'coriolisterm':'c',
                     'gravityterm':'g',
                     'frictionterm':'f',
                     'baseregressor':'Hb'}
    for i, o in enumerate(outputs):
        if o in outnames_dict:
            outputs[i] = outnames_dict[o]
    outputs_set = set(outputs)

    se = Subexprs()
    exprs = {}

    fprint('Computing geometric model ... ')
    rbtgeo = Geometry(rbtdef, se.collect)
    fprint('done\n')

    if 'T' in outputs_set:
        exprs['T'] = rbtgeo.T[-1]

    if 'J' in outputs_set:
        fprint('Computing kinematic model ... ')
        rbtkin = Kinematics(rbtdef, rbtgeo, se.collect)
        fprint('done\n')
        exprs['J'] = rbtkin.J[-1]

    rbtdyn = Dynamics(rbtdef, rbtgeo)

    if 'tau' in outputs_set:
        fprint('Generating inverse dynamics ... ')
        rbtdyn.gen_invdyn(se.collect)
        exprs['tau'] = rbtdyn.invdyn
        fprint('done\n')

    if 'M' in outputs_set:
        fprint('Generating inertia matrix ... ')
        rbtdyn.gen_inertiamatrix(se.collect)
        exprs['M'] = rbtdyn.M
        fprint('done\n')

    if 'c' in outputs_set:
        fprint('Computing coriolis term ... ')
        rbtdyn.gen_coriolisterm(se.collect)
        exprs['c'] = rbtdyn.c
        fprint('done\n')

    if 'g' in outputs_set:
        fprint('Computing gravity term ... ')
        rbtdyn.gen_gravityterm(se.collect)
        exprs['g'] = rbtdyn.g
        fprint('done\n')

    if 'f' in outputs_set:
        fprint('Computing friction term ... ')
        rbtdyn.gen_frictionterm(se.collect)
        exprs['f'] = rbtdyn.f
        fprint('done\n')

    if {'H', 'Hb'} & outputs_set or 'base dynamic parameters' in defs_dict['generate']:
        fprint('Computing regressor ... ')
        rbtdyn.gen_regressor(se.collect)
        exprs['H'] = rbtdyn.H
        fprint('done\n')

    if 'Hb' in outputs_set or 'base dynamic parameters' in defs_dict['generate']:

        fprint('Generating regressor function ... ')
        func_def_regressor = robot_code_to_func('python', se.get(exprs['H']),
                                                'H', 'regressor_func', rbtdef)
        # global sin, cos, sign
        # sin = numpy.sin
        # cos = numpy.cos
        # sign = numpy.sign
        l = locals()
        exec_(func_def_regressor, globals(), l)
        regressor_func = l['regressor_func']
        fprint('done\n')

        fprint('Computing base parameters ... ')
        rbtdyn.calc_base_parms(regressor_func)
        fprint('done\n')

        if 'Hb' in outputs_set:
            fprint('Computing base regressor ... ')
            exprs['Hb'] = exprs['H'] * rbtdyn.Pb
            fprint('done\n')

    fprint('Generating source code file ... ')
    lang = source['lang'].lower()
    lang = {'python':'py', 'c++':'c'}.get(lang, lang)
    code = se.get([exprs[o] for o in outputs])
    srccode = robot_code_to_func(lang, code, outputs, source['funcname'], rbtdef)
    if outputs_set & {'tau', 'M', 'g', 'c', 'f', 'H'}:
        comment = {'py':'# ', 'c':'// '}.get(lang, '')
        parms_str =  ", ".join([str(p) for p in rbtdef.dynparms()])
        parms_str = textwrap.wrap(parms_str, 77)
        parms_str = comment + ('\n' + comment).join(parms_str)
        parms_str = '\n' + comment + ' --- Order of dynamic parameters ---\n' \
            + parms_str
        parms_str += '\n\n'
        srccode = parms_str + srccode
    fprint('done\n')

    fprint("Saving source code into '%s' file ... "%source['filename'])
    f = open(source['filename'], 'w')
    f.write(srccode)
    f.close()
    fprint('done\n')

    if 'base dynamic parameters' in defs_dict['generate']:
        filename = defs_dict['generate']['base dynamic parameters']['filename']
        fprint("Saving base parameter info into '%s' file ... "%filename)
        baseparms = 'Base parameters:\n'
        for i, b in enumerate(rbtdyn.baseparms):
            baseparms += 'b%d = %s\n'%(i+1, b.n())

        def print_mat(m, f='%g'):
            s = ''
            for i in range(m.rows):
                s += ' '.join([f%m[i,j] for j in range(m.cols)])
                s += '\n'
            return s

        baseparms += '\nPb:\n'
        baseparms += print_mat(rbtdyn.Pb)
        baseparms += '\nPd:\n'
        baseparms += print_mat(rbtdyn.Pd)
        baseparms += '\nKd:\n'
        baseparms += print_mat(rbtdyn.Kd, '%.10g')

        baseparms += '\n'
        f = open(filename, 'w')
        f.write(baseparms)
        f.close()
        fprint('done\n')

    return



example='''%YAML 1.2
---
robot:
    name: 3-DOF Planar Robot

    DH parameters:
        - [0, 3.0, 0, q]                   # < order: alpha, a, d, theta
        - [0, 2, 0, q + pi]                #      q means the joint position
        - [0, 2.5, 0, q + pi/2]

    DH convention: standard                # < standard OR modified

    extra:                                 # < OPTIONAL definitions
        friction model: simple             # < simple OR None (None by default)
        gravity acceleration: [0, 0, -9.8] # < in base frame coords
                                           #     ([0, 0, -9.8] by default)

generate:
    source:
        filename: planar3.c
        lang: c
        funcname: rbt
        outputs: [T, M, g, c]   # < options: T, J, tau, M, c, g, f, H, Hb
                                #     corresponding, respectively, to:
                                #      end-effector transformation, Jacobian,
                                #   optional   inverse dynamics (no friction), inertia
                                #      matrix, Coriolis term, gravity term,
                                #      friction term, regressor, base regressor

    base dynamic parameters:    # < OPTIONAL
        filename: planar3.txt   # <
...
'''




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate robot model from a given .yml definition file.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('filename', nargs='?', help='robot definitions YAML file')
    group.add_argument('--example', action='store_true', help='print example file')
    args = parser.parse_args()

    if args.example:
        print(example)
    else:
        f = open(args.filename)
        d = yaml.safe_load(f)
        f.close()
        gen_robot(d)
