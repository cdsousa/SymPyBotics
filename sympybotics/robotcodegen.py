
import sympy

from . import symcode


def _gen_q_dq_ddq_subs(rbtdef):
    subs = {}
    for i in reversed(range(rbtdef.dof)):
        subs[rbtdef.ddq[i]] = 'ddq[' + str(i) + ']'
        subs[rbtdef.dq[i]] = 'dq[' + str(i) + ']'
        subs[rbtdef.q[i]] = 'q[' + str(i) + ']'
    return subs


def _gen_parms_subs(parms_symbols, name='parms'):
    subs = {}
    for i, parm in enumerate(parms_symbols):
        subs[parm] = name + '[' + str(i) + ']'
    return subs


def dyn_code_to_func(lang, code, funcname, rbtdef):
    func_parms = []
    subs_pairs = {}

    code_symbols = set()
    for iv, se in code[0]:
        code_symbols.update(se.free_symbols)
    for e in code[1]:
        code_symbols.update(e.free_symbols)

    if not code_symbols.isdisjoint(set(rbtdef.dynparms())):
        func_parms.append('parms')
        subs_pairs.update(_gen_parms_subs(rbtdef.dynparms(), 'parms'))

    qderivlevel = -1
    if not code_symbols.isdisjoint(set(rbtdef.q)):
        qderivlevel = 0
    if not code_symbols.isdisjoint(set(rbtdef.dq)):
        qderivlevel = 1
    if not code_symbols.isdisjoint(set(rbtdef.ddq)):
        qderivlevel = 2

    if qderivlevel >= 0:
        for i in range(qderivlevel + 1):
            func_parms.append('d'*i + 'q')
        subs_pairs.update(_gen_q_dq_ddq_subs(rbtdef))

    return symcode.generation.code_to_func(lang, code, funcname, func_parms,
                                           subs_pairs)
