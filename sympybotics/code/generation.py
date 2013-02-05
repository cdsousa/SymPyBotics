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



import copy
import string
import sympy
import sys
import re


options = {}
options['unroll_square'] = True


def code_back_to_exprs(code):
    
    ivars = copy.deepcopy(code[0])
    exps = copy.deepcopy(code[1])
    
    for i in range(len(ivars)):
        
        Dxreplace = { ivars[i][0] : ivars[i][1] }
        
        for j in range(i+1,len(ivars)):
                if ivars[j][1].has(ivars[i][0]):
                    ivars[j] = ( ivars[j][0], ivars[j][1].xreplace( Dxreplace ) )
        
        for j in range(len(exps)):
                if sympy.sympify(exps[j]).has(ivars[i][0]):
                    exps[j] = exps[j].xreplace( Dxreplace )
    
    return exps

def _ccode( expr, ):
  code = sympy.ccode( expr )
  if options['unroll_square']:
    return re.sub(r'pow\(([^,]*), 2\)', r'((\1)*(\1))', code)
  else:
    return code
  
def code_to_string( code, outvar_name='out', indent='', realtype='', line_end='' ):
    
    codestr = ''
    
    if realtype: realtype += ' '
    
    for i in range( len(code[0]) ) :
        codestr += indent + realtype + sympy.ccode( code[0][i][0] ) + ' = ' + _ccode( code[0][i][1] ) + line_end + '\n'
    
    codestr += '\n'
    for i in range( len(code[1]) ) :
        codestr += indent + outvar_name + '['+str(i)+'] = ' + _ccode( code[1][i] ) + line_end + '\n'
    
    return codestr

def codestring_count( codestring, resume=False ):
  ops = []
  ops += [( '=' , int(codestring.count('=')) )]
  ops += [( '+' , int(codestring.count('+')) )]
  ops += [( '-' , int(codestring.count('-')) )]
  ops += [( '*' , int(codestring.count('*')) )]
  ops += [( '/' , int(codestring.count('/')) )]
  ops += [( 'pow' , int(codestring.count('pow')) )]
  ops += [( 'sin' , int(codestring.count('sin')) )]
  ops += [( 'cos' , int(codestring.count('cos')) )]
  if not resume:
    return ops
  else:
    adds = int(codestring.count('+'))+int(codestring.count('-'))
    muls = int(codestring.count('*'))+int(codestring.count('/'))+int(codestring.count('pow'))
    return ops, {'add':adds, 'mul':muls, 'total':adds+muls }


def gen_py_func( code, func_parms, subs_pairs, func_name='func', outvar_name='out' ):

    indent = 4*' '

    pycode = 'def ' + func_name + '('
    if func_parms:
        pycode += ' ' + func_parms[0]
        for parm in func_parms[1:] :
            pycode += ', ' + parm
        pycode += ' '
    pycode += ') :\n\n'

    pycode += indent + outvar_name + ' = [0]*' + str( len(code[1]) ) + '\n\n'

    mainpycode = code_to_string( code, outvar_name, indent )
    for (old,new) in subs_pairs: mainpycode = mainpycode.replace(old,new)
    pycode += mainpycode

    pycode += '\n' + indent + 'return ' + outvar_name + '\n'
    return pycode

def gen_c_func( code, func_parms, subs_pairs, func_name='func', outvar_name='out' ):

    indent = 2*' '

    ccode = 'void ' + func_name + '( double* ' + outvar_name
    for parm in func_parms :
        ccode += ', const double* ' + parm
    ccode += ' )\n{\n'

    mainccode = code_to_string( code, outvar_name, indent, 'double', ';' )
    for (old,new) in subs_pairs: mainccode = mainccode.replace(old,new)

    ccode += mainccode + '\n'+indent+'return;\n}\n'
    return ccode

def gen_pyx_func( code, func_parms, subs_pairs, func_name='func', outvar_name='out' ):

    indent = 4*' '

    ccode = 'cdef void ' + func_name + '( double* ' + outvar_name
    for parm in func_parms :
        ccode += ', double* ' + parm
    ccode += ' ):\n'

    mainccode = code_to_string( code, outvar_name, indent, 'cdef double' )
    for (old,new) in subs_pairs: mainccode = mainccode.replace(old,new)

    ccode += mainccode + '\n'+indent+'return\n'
    return ccode


def code_to_func( lang, code, func_name, func_parms, subs_pairs ):
  lang = lang.lower()
  if lang in ['python','py'] : gen_func = gen_py_func
  elif lang in ['cython','pyx'] : gen_func = gen_pyx_func
  elif lang in ['c','c++'] : gen_func = gen_c_func
  else: raise Exception('chosen language not supported.')
  return gen_func( code, func_parms, subs_pairs, func_name, func_name+'_out' )


