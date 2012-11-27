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
import random
import string
import sympy


def code_apply_func( code, func ):
  retcode = copy.deepcopy(code)
  for i in range(len(retcode[0])):
    retcode[0][i] = ( retcode[0][i][0], func(retcode[0][i][1]) )
  for i in range(len(retcode[1])):
    retcode[1][i] = func(retcode[1][i])
  return retcode


def code_subs( code, subs_dict ):
  return code_apply_func( code, lambda x: x.subs( subs_dict ) )

def code_cse( code, auxvarname = None ):

    if not auxvarname:
      auxvarname = 'cse_' + ''.join([random.choice(string.ascii_lowercase+string.digits) for _ in range(4)]) + '_ivar_'
      
    code_in = code

    codecse = sympy.cse( sympy.sympify( [i[1] for i in code_in[0]] + code_in[1] ), \
                         sympy.cse_main.numbered_symbols(auxvarname) )
    
    auxv1_num = len(code_in[0])
    
    if auxv1_num == 0:
        return codecse
    
    A1 = list(zip( list(zip(*code_in[0]))[0] , codecse[1][:auxv1_num] ))
    Acse = codecse[0]
    
    def getdeps(expr,assignments):
        out = []
        search_assignments = copy.copy(assignments)
        for assign in search_assignments:
            if sympy.sympify(expr).has(assign[0]) and assign in assignments:
                assignments.remove(assign)
                out += getdeps(assign[1],assignments) + [assign]
        return out
    
    codemerge = []
    for a1 in A1:
        codemerge += getdeps(a1[1],Acse)
        codemerge.append( a1 )
    for acse in Acse:
        codemerge.append( acse )
    
    retcode = ( codemerge , codecse[1][auxv1_num:] )
    
    return retcode


def code_remove_not_compound( code ):
    
    debug = False
    removed=0
    
    retcode = copy.deepcopy(code)
    toremove = []
    
    for i in range(len(retcode[0])):
        v = retcode[0][i][0]
        e = retcode[0][i][1]
        if e.is_Atom or (-e).is_Atom:
            if debug: print(i,v,e,'is atom')
            retcode = code_subs( retcode, {v:e} )
            toremove.append(i)
            removed += 1
            if debug: print('  poped')

    for r in reversed(toremove):
        retcode[0].pop(r)
    
    if debug: print('removed',removed)
    return retcode


def code_remove_not_or_once_used( code ):
    
    debug = False
    removed=0
    
    retcode = copy.deepcopy(code)
    
    for i in range(len(retcode[0])-1,-1,-1):
        v = retcode[0][i][0]
        e = retcode[0][i][1]
        uses = 0
        usedin_ai = None
        usedin_oi = None
        if debug: print(i,v)
        for ai in range(i+1,len(retcode[0])):
            count = sympy.sympify(retcode[0][ai][1]).has(v)
            if count:
                uses += count
                usedin_ai = ai
                # if debug: print('  used in',ai)
        for oi in range(len(retcode[1])):
            count = sympy.sympify(retcode[1][oi]).count(v)
            if count:
                uses += count
                usedin_oi = oi
                # if debug: print('  used in out',oi)
        
        if uses == 0:
            retcode[0].pop(i)
            removed += 1
            if debug: print('  not used: removed')
        elif uses == 1:
            if usedin_ai != None:
                retcode[0][usedin_ai] = ( retcode[0][usedin_ai][0] , retcode[0][usedin_ai][1].subs({v:e}) )
                if debug: print('  used once: removed and substituted in',retcode[0][usedin_ai][0])
            elif usedin_oi != None:
                retcode[1][usedin_oi] = retcode[1][usedin_oi].subs({v:e})
                if debug: print('  used once: removed and substituted in out',usedin_oi)
            retcode[0].pop(i)
            removed += 1
    
    if debug: print('removed',removed)
    return retcode


def code_rename_ivars_unsafe(code, ivarnames ):
    
  retcode = copy.deepcopy(code)
   
  for i in range(len(retcode[0])):

    new_symbol = sympy.Symbol(ivarnames+str(i),real=True)
    Dsubs = { retcode[0][i][0] : new_symbol }
  
    retcode[0][i] = ( new_symbol , retcode[0][i][1] )
        
    for j in range(i+1,len(retcode[0])):
      retcode[0][j] = ( retcode[0][j][0], retcode[0][j][1].subs( Dsubs ) )
  
    for j in range(len(retcode[1])):
      retcode[1][j] = sympy.sympify(retcode[1][j]).subs( Dsubs )
  
  return retcode


def code_make_output_single_vars(code, ivarnames=None ):

    retcode = copy.deepcopy(code)

    if ivarnames:
        cnt = 0

    else:
        if retcode[0]:
            lastivar = str(retcode[0][-1][0])
            for i in range(len(lastivar)):
                if not lastivar[-i-1].isdigit(): break
            if i > 0:
                cnt = int(lastivar[-i:]) + 1
                ivarnames = lastivar[:-i]
            else:
                cnt = 1
                ivarnames = lastivar + '_'
        else:
            cnt = 0
            ivarnames = 'outputiv_'

    for i in range(len(retcode[1])):
        if not sympy.sympify(retcode[1][i]).is_Atom:
            new_symbol = sympy.Symbol(ivarnames+str(cnt),real=True)
            retcode[0].append( (new_symbol, retcode[1][i]) )
            retcode[1][i] = new_symbol
            cnt += 1

    return retcode


def optimize_code( code, ivarnames='iv', singlevarout=False, debug = True ) :
  if debug: print('Optimizing code')
  retcode = copy.deepcopy(code)
  if debug: print('code_apply_func trigsimp')
  retcode = code_apply_func( retcode, lambda x: sympy.trigsimp(x) )
  if debug: print('code_remove_not_compound')
  retcode = code_remove_not_compound(retcode)
  if debug: print('code_remove_not_or_once_used')
  retcode = code_remove_not_or_once_used(retcode)
  if debug: print('code_cse')
  retcode = code_cse(retcode,'cse')
  #if debug: print('code_remove_not_compound')
  #retcode = code_remove_not_compound(retcode)
  #if debug: print('code_remove_not_or_once_used')
  #retcode = code_remove_not_or_once_used(retcode)
  #if debug: print('code_cse 2')
  #retcode = code_cse(retcode,'cse2')
  if debug: print('code_remove_not_compound')
  retcode = code_remove_not_compound(retcode)
  if debug: print('code_remove_not_or_once_used')
  retcode = code_remove_not_or_once_used(retcode)
  if debug: print('code_rename_ivars (unsafe)')
  retcode = code_rename_ivars_unsafe(retcode, ivarnames=ivarnames)
  if singlevarout:
    if debug: print('code_make_output_single_vars')
    retcode = code_make_output_single_vars(retcode)
  if debug: print('Done.')
  return retcode


def code_back_to_expressions(code):
    
    ivars = copy.deepcopy(code[0])
    exps = copy.deepcopy(code[1])
    
    for i in range(len(ivars)):
        
        Dsubs = { ivars[i][0] : ivars[i][1] }
        
        for j in range(i+1,len(ivars)):
                if ivars[j][1].has(ivars[i][0]):
                    ivars[j] = ( ivars[j][0], ivars[j][1].subs( Dsubs ) )
        
        for j in range(len(exps)):
                if sympy.sympify(exps[j]).has(ivars[i][0]):
                    exps[j] = exps[j].subs( Dsubs )
    
    return exps


def code_to_string( code, outvar_name='out', indent='', realtype='', line_end='' ):
    
    codestr = ''
    
    if realtype: realtype += ' '
    
    for i in range( len(code[0]) ) :
        codestr += indent + realtype + sympy.ccode( code[0][i][0] ) + ' = ' + sympy.ccode( code[0][i][1] ) + line_end + '\n'
    
    codestr += '\n'
    for i in range( len(code[1]) ) :
        codestr += indent + outvar_name + '['+str(i)+'] = ' + sympy.ccode( code[1][i] ) + line_end + '\n'
    
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


def sympymatrix_to_func( lang, ivars, matrix, func_name, func_parms, subs_pairs ):
    code = optimize_code( (ivars, sympy.flatten(matrix)), ivarnames='aux' )
    return code_to_func( lang, code, func_name, func_parms, subs_pairs )

