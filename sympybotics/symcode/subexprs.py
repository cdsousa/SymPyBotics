import sympy
import sympy.utilities
from sympy.simplify.cse_main import cse_optimizations, preprocess_for_cse, \
    postprocess_for_cse

import collections


class Subexprs(object):

    def __init__(self, optimizations=None, postprocess=None):

        if optimizations is None:
            # Pull out the default here just in case there are some weird
            # manipulations of the module-level list in some other thread.
            optimizations = list(cse_optimizations)
        self._optimizations = optimizations
        self._postprocess = postprocess

        self._tmp_symbols = sympy.utilities.iterables.numbered_symbols(
            'tmp', start=0, real=True)

        self._subexp_iv = dict()
        self._commutatives = dict()

    class _ordered_len(object):

        def __init__(self):
            self.lenidxs = [0]

        def insert(self, list, item):
            l = len(item)
            over = l - (len(self.lenidxs) - 1)
            if over > 0:
                self.lenidxs.extend([self.lenidxs[-1]] * over)
            list.insert(self.lenidxs[l], item)
            for j in range(l, len(self.lenidxs)):
                self.lenidxs[j] += 1

        def pop(self, list, index):
            l = len(list.pop(index))
            for j in range(l, len(self.lenidxs)):
                self.lenidxs[j] -= 1

    def _parse_commutative(self, expr):

        exprtype = type(expr)
        args_input = set(expr.args)

        if exprtype not in self._commutatives:
            argsets = []
            argset_orderlens = self._ordered_len()
            argset_orderlens.insert(argsets, args_input)
            self._commutatives[exprtype] = (argsets, argset_orderlens)
            ivar = next(self._tmp_symbols)
            self._subexp_iv[expr] = ivar
            return ivar

        argsets = self._commutatives[exprtype][0]
        argset_orderlens = self._commutatives[exprtype][1]

        ivar = None
        args_to_remove = []
        args_to_insert = []

        # for 2 args input exprs, bypass comparison with other 2 args exprs
        if len(args_input) == 2:
            init = argset_orderlens.lenidxs[2]
        else:
            init = 0

        for i in range(init, len(argsets)):
            args_other = argsets[i]

            com = args_input.intersection(args_other)
            if len(com) > 1:

                diff_args_input = args_input.difference(com)
                diff_args_other = args_other.difference(com)

                if not diff_args_input:
                    # args_input is strict subset of args_other

                    ivar = next(self._tmp_symbols)
                    self._subexp_iv[exprtype(*args_input)] = ivar
                    args_to_insert.append(args_input)

                    args_other = diff_args_other
                    args_other.add(ivar)
                    self._subexp_iv[exprtype(*args_other)] = \
                        self._subexp_iv.pop(exprtype(*argsets[i]))
                    args_to_remove.append(i)
                    args_to_insert.append(args_other)

                    break

                elif not diff_args_other:
                    # args_other is strict subset of args_input

                    args_input = diff_args_input
                    args_input.add(self._subexp_iv[exprtype(*args_other)])

                    ivar = self._subexp_iv.get(exprtype(*args_input), None)

                    if ivar or len(args_input) == 2:
                        break

                else:  # args_input != com != args_other

                    ivar_com = next(self._tmp_symbols)
                    self._subexp_iv[exprtype(*com)] = ivar_com
                    args_to_insert.append(com)  # argsets.append(com)

                    args_other = diff_args_other
                    args_other.add(ivar_com)
                    self._subexp_iv[exprtype(*args_other)] = \
                        self._subexp_iv.pop(exprtype(*argsets[i]))
                    args_to_remove.append(i)
                    args_to_insert.append(args_other)

                    args_input = diff_args_input
                    args_input.add(ivar_com)

                    if len(args_input) == 2:
                        break

        if ivar is None:
            ivar = next(self._tmp_symbols)
            self._subexp_iv[exprtype(*args_input)] = ivar
            args_to_insert.append(args_input)

        for i in reversed(sorted(args_to_remove)):
            argset_orderlens.pop(argsets, i)
        for args in args_to_insert:
            argset_orderlens.insert(argsets, args)

        return ivar

    def _parse(self, expr):

        if expr.is_Atom:
            # Exclude atoms, since there is no point in renaming them.
            return expr

        if sympy.iterables.iterable(expr):
            return expr

        subexpr = type(expr)(*map(self._parse, expr.args))

        if subexpr in self._subexp_iv:
            return self._subexp_iv[subexpr]

        if subexpr.is_Mul or subexpr.is_Add:
            return self._parse_commutative(subexpr)
        else:
            ivar = next(self._tmp_symbols)
            self._subexp_iv[subexpr] = ivar
            return ivar

    def collect(self, exprs):

        if isinstance(exprs, sympy.Basic):  # if only one expression is passed
            exprs = [exprs]
            is_single_expr = True
        else:
            is_single_expr = False

        # Preprocess the expressions to give us better optimization
        # opportunities.
        prep_exprs = [preprocess_for_cse(e, self._optimizations)
                      for e in exprs]

        out_exprs = map(self._parse, prep_exprs)

        if is_single_expr:
            return out_exprs[0]
        elif isinstance(exprs, sympy.Matrix):
            return sympy.Matrix(exprs.rows, exprs.cols, out_exprs)
        else:
            return out_exprs

    def get(self, exprs=None, symbols=None):

        if symbols is None:
            symbols = sympy.utilities.iterables.numbered_symbols()
        else:
            # In case we get passed an iterable with an __iter__ method
            # instead of an actual iterator.
            symbols = iter(symbols)

        if isinstance(exprs, sympy.Basic):  # if only one expression is passed
            exprs = [exprs]

        # Find all of the repeated subexpressions.

        ivar_se = {iv: se for se, iv in self._subexp_iv.iteritems()}

        used_ivs = set()
        repeated = set()

        def _find_repeated_subexprs(subexpr):
            if subexpr.is_Atom:
                symbs = [subexpr]
            else:
                symbs = subexpr.args
            for symb in symbs:
                if symb in ivar_se:
                    if symb not in used_ivs:
                        _find_repeated_subexprs(ivar_se[symb])
                        used_ivs.add(symb)
                    else:
                        repeated.add(symb)

        for expr in exprs:
            _find_repeated_subexprs(expr)

        # Substitute symbols for all of the repeated subexpressions.
        # remove temporary replacements that weren't used more than once

        tmpivs_ivs = dict()
        ordered_iv_se = collections.OrderedDict()

        def _get_subexprs(args):
            args = list(args)
            for i, symb in enumerate(args):
                if symb in ivar_se:
                    if symb in tmpivs_ivs:
                        args[i] = tmpivs_ivs[symb]
                    else:
                        subexpr = ivar_se[symb]
                        subexpr = type(subexpr)(*_get_subexprs(subexpr.args))
                        if symb in repeated:
                            ivar = next(symbols)
                            ordered_iv_se[ivar] = subexpr
                            tmpivs_ivs[symb] = ivar
                            args[i] = ivar
                        else:
                            args[i] = subexpr
            return args

        out_exprs = _get_subexprs(exprs)

        # Postprocess the expressions to return the expressions to canonical
        # form.
        ordered_iv_se_notopt = ordered_iv_se
        ordered_iv_se = collections.OrderedDict()
        for i, (ivar, subexpr) in enumerate(ordered_iv_se_notopt.items()):
            subexpr = postprocess_for_cse(subexpr, self._optimizations)
            ordered_iv_se[ivar] = subexpr
        out_exprs = [postprocess_for_cse(e, self._optimizations)
                     for e in out_exprs]

        if isinstance(exprs, sympy.Matrix):
            out_exprs = sympy.Matrix(exprs.rows, exprs.cols, out_exprs)
        if self._postprocess is None:
            return ordered_iv_se.items(), out_exprs
        return self._postprocess(ordered_iv_se.items(), out_exprs)


def fast_cse(exprs, symbols='aux'):
    se = Subexprs()
    return se.get(se.collect(exprs))


class WholeSubexprs(object):

    def __init__(self, *args, **kwargs):

        self._tmp_symbols = sympy.utilities.iterables.numbered_symbols()
        self._subexp_iv = list()

    def collect(self, exprs):

        if isinstance(exprs, sympy.Basic):  # if only one expression is passed
            exprs = [exprs]
            is_single_expr = True
        else:
            is_single_expr = False

        out_exprs = []
        for expr in exprs:
            if expr.is_Atom or (-expr).is_Atom:
                out_exprs.append(expr)
            else:
                iv = next(self._tmp_symbols)
                self._subexp_iv.append((iv, expr))
                out_exprs.append(iv)

        if is_single_expr:
            return out_exprs[0]
        elif isinstance(exprs, sympy.Matrix):
            return sympy.Matrix(exprs.rows, exprs.cols, out_exprs)
        else:
            return out_exprs

    def get(self, exprs=None, symbols=None):
        return self._subexp_iv, exprs
