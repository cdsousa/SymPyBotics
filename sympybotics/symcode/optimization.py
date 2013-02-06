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
import sympy
import sys



def apply_func( code, func ):
  sube = [None]*len(code[0])
  oute = copy.deepcopy(code[1])
  for (i,(v,e)) in enumerate(code[0]):
    sube[i] = ( v, func(e) )
  for i,e in enumerate(code[1]):
    oute[i] = func(e)
  return (sube, oute)


def xreplace( code, xreplace_dict ):
  return apply_func( code, lambda x: x.xreplace( xreplace_dict ) )


def optim_dce( code ):
    """Performe 'dead code elimination' optimization on code."""
    
    ivs = []
    
    n = len(code[0])
    
    for i in range(len(code[0])-1,-1,-1):
        s = code[0][i][0]
        
        used = False
        
        for _,e in ivs:
            if sympy.sympify(e).has(s):
              used = True
              break
        if not used:
            for e in code[1]:
              if sympy.sympify(e).has(s):
                used = True
                break
        
        if used:
          ivs.append(code[0][i])
          
    return list(reversed(ivs)), code[1][:]


def optim_cse( code, auxvarname = 'cse' ):
    """Performe 'common sub-expressions elimination' optimization on code."""

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


def optim_cp( code ):
    """Performe 'copy propagation' optimization on code."""
    
    debug = False
    removed=0
    
    retcode = copy.deepcopy(code)
    toremove = []
    
    for i in range(len(retcode[0])):
        v = retcode[0][i][0]
        e = retcode[0][i][1]
        if e.is_Atom:
            if debug: print(i,v,e,'is atom')
            retcode = xreplace( retcode, {v:e} )
            toremove.append(i)
            removed += 1
            if debug: print('  poped')

    for r in reversed(toremove):
        retcode[0].pop(r)
    
    if debug: print('removed',removed)
    return retcode


def optim_dce_sup( code ):
    """Performe 'dead code elimination' and 'single use propagation' optimizations on code."""
    
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
                retcode[0][usedin_ai] = ( retcode[0][usedin_ai][0] , retcode[0][usedin_ai][1].xreplace({v:e}) )
                if debug: print('  used once: removed and substituted in',retcode[0][usedin_ai][0])
            elif usedin_oi != None:
                retcode[1][usedin_oi] = retcode[1][usedin_oi].xreplace({v:e})
                if debug: print('  used once: removed and substituted in out',usedin_oi)
            retcode[0].pop(i)
            removed += 1
    
    if debug: print('removed',removed)
    return retcode


def rename_ivars_unsafe(code, ivarnames ):
    
  retcode = copy.deepcopy(code)
   
  for i in range(len(retcode[0])):

    new_symbol = sympy.Symbol(ivarnames+str(i),real=True)
    Dxreplace = { retcode[0][i][0] : new_symbol }
  
    retcode[0][i] = ( new_symbol , retcode[0][i][1] )
        
    for j in range(i+1,len(retcode[0])):
      retcode[0][j] = ( retcode[0][j][0], retcode[0][j][1].xreplace( Dxreplace ) )
  
    for j in range(len(retcode[1])):
      retcode[1][j] = sympy.sympify(retcode[1][j]).xreplace( Dxreplace )
  
  return retcode


def make_output_single_vars(code, ivarnames=None ):

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

def _fprint(x):
  print(x)
  sys.stdout.flush()

def fully_optimize_code( code, ivarnames=None, singlevarout=False, clearcache=0, debug = True ) :
  
  if debug: _fprint('Optimizing code')
  

  if debug: _fprint('dead code elimination and single use propagation')
  code = optim_dce_sup(code)
  if clearcache > 1: sympy.cache.clear_cache()
  
  if debug: _fprint('common sub-expressions elimination')
  code = optim_cse(code,'cse')
  if clearcache > 1: sympy.cache.clear_cache()
  
  if debug: _fprint('copy propagation')
  code = optim_cp(code)
  if clearcache > 1: sympy.cache.clear_cache()
  
  if debug: _fprint('dead code elimination and single use propagation')
  code = optim_dce_sup(code)
  if clearcache > 1: sympy.cache.clear_cache()
  
  if ivarnames:
    if debug: _fprint('code_rename_ivars (unsafe)')
    code = rename_ivars_unsafe(code, ivarnames=ivarnames)
    if clearcache > 1: sympy.cache.clear_cache()
  
  if singlevarout:
    if debug: _fprint('code_make_output_single_vars')
    code = make_output_single_vars(code)
    
  if clearcache: sympy.cache.clear_cache()
  
  if debug: print('Done.')
  
  return code


