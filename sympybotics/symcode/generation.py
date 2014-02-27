
import copy
import sympy
import re


options = {}
options['unroll_square'] = True


def apply_func(code, func, apply_to_ivs=True):
    if apply_to_ivs:
        code_ivs = [(func(iv), func(se)) for iv, se in code[0]]
    else:
        code_ivs = [(iv, func(se)) for iv, se in code[0]]
    code_exprs = []
    for expr in code[1]:
        if isinstance(expr, sympy.MatrixBase):
            expr = expr.applyfunc(func)
        else:
            expr = func(expr)
        code_exprs.append(expr)
    return code_ivs, code_exprs


def xreplace(code, xreplace_dict):
    return apply_func(code, lambda x: x.xreplace(xreplace_dict))


def code_back_to_exprs(code):

    ivars = copy.deepcopy(code[0])
    exps = copy.deepcopy(code[1])

    for i in range(len(ivars)):

        Dxreplace = {ivars[i][0]: ivars[i][1]}

        for j in range(i + 1, len(ivars)):
            if ivars[j][1].has(ivars[i][0]):
                ivars[j] = (ivars[j][0], ivars[j][1].xreplace(Dxreplace))

        for j in range(len(exps)):
            if sympy.sympify(exps[j]).has(ivars[i][0]):
                exps[j] = exps[j].xreplace(Dxreplace)

    return exps


def _ccode(expr, ):
    code = sympy.ccode(expr)
    if options['unroll_square']:
        return re.sub(r'pow\(([^,]*), 2\)', r'((\1)*(\1))', code)
    else:
        return code


def _juliacode(expr, ):
    code = sympy.printing.lambdarepr.lambdarepr(expr)
    return code.replace('**', '^')


def code_to_string(code, out_parms, printer, indent='', realtype='',
                   line_end='', outidxoffset=0):

    codestr = ''

    if realtype:
        realtype += ' '

    for i in range(len(code[0])):
        codestr += indent + realtype + \
            printer(code[0][i][0]) + ' = ' + printer(
                code[0][i][1]) + line_end + '\n'

    for c, out in enumerate(out_parms):
        codestr += '\n'
        for i in range(len(code[1][c])):
            codestr += indent + out + \
                '[' + str(i+outidxoffset) + '] = ' + \
                printer(code[1][c][i]) + line_end + '\n'

    return codestr


def codestring_count(codestring, resume=False):
    ops = []
    ops += [('=', int(codestring.count('=')))]
    ops += [('+', int(codestring.count('+')))]
    ops += [('-', int(codestring.count('-')))]
    ops += [('*', int(codestring.count('*')))]
    ops += [('/', int(codestring.count('/')))]
    ops += [('pow', int(codestring.count('pow')))]
    ops += [('sin', int(codestring.count('sin')))]
    ops += [('cos', int(codestring.count('cos')))]
    if not resume:
        return ops
    else:
        adds = int(codestring.count('+')) + int(codestring.count('-'))
        muls = int(codestring.count('*')) + int(
            codestring.count('/')) + int(codestring.count('pow'))
        return ops, {'add': adds, 'mul': muls, 'total': adds + muls}


def gen_py_func(code, out_parms, func_parms, func_name='func'):

    indent = 4 * ' '

    pycode = 'def ' + func_name + '('
    pycode += ', '.join(func_parms)
    pycode += '):\n\n'

    for i, out in enumerate(out_parms):
        pycode += indent + out + ' = [0]*' + str(len(code[1][i])) + '\n'

    pycode += '\n'

    mainpycode = code_to_string(code, out_parms,
                                sympy.printing.lambdarepr.lambdarepr, indent)

    pycode += mainpycode

    pycode += '\n' + indent + 'return '
    pycode += ', '.join(out_parms)

    pycode = pycode.replace('\n\n', '\n#\n')

    return pycode


def gen_c_func(code, out_parms, func_parms, func_name='func'):

    indent = 2 * ' '

    ccode = 'void ' + func_name + '( double* '

    ccode += ', double* '.join(out_parms)

    ccode += ', const double* '
    ccode += ', const double* '.join(func_parms)

    ccode += ' )\n{\n'

    mainccode = code_to_string(code, out_parms, _ccode, indent, 'double', ';')

    ccode += mainccode + '\n' + indent + 'return;\n}'

    ccode = ccode.replace('\n\n', '\n//\n')

    return ccode


def gen_julia_func(code, out_parms, func_parms, func_name='func'):

    indent = 4 * ' '

    ccode = 'function ' + func_name + '('

    ccode += ', '.join(out_parms)

    ccode += ', '
    ccode += ', '.join(func_parms)

    ccode += ')\n\n'

    code = code[0], [(e.T if isinstance(e, sympy.MatrixBase) else e)
                     for e in code[1]]

    mainccode = code_to_string(code, out_parms, _juliacode, indent, '', '',
                               outidxoffset=1)

    ## pass from 0-idexed to 1-indexed arrays
    #mainccode = re.sub(
        #r"\[([0-9]+)\]",
        #lambda m: '[' + str(int(m.group(1))+1) + ']',
        #mainccode)

    ccode += mainccode + '\n' + indent + 'return '
    ccode += ', '.join(out_parms)
    ccode += '\nend'

    ccode = ccode.replace('\n\n', '\n#\n')

    return ccode


def code_to_func(lang, code, out_parms, func_name, func_parms, symb_replace):

    if not isinstance(code[1], list):
        code = (code[0], [code[1]])
    if not isinstance(out_parms, list):
        out_parms = [out_parms]

    lang = lang.lower()
    if lang in ['python', 'py']:
        gen_func = gen_py_func
    elif lang in ['c', 'c++']:
        gen_func = gen_c_func
    elif lang in ['julia', 'jl']:
        gen_func = gen_julia_func
    else:
        raise Exception('chosen language not supported.')

    if symb_replace:
        sympified_replace = {}
        for k, v in symb_replace.items():
            if isinstance(k, str):
                k = sympy.Symbol(k)
            if isinstance(v, str):
                v = sympy.Symbol(v)
            sympified_replace[k] = v
        code = xreplace(code, sympified_replace)

    return gen_func(code, out_parms, func_parms, func_name)
