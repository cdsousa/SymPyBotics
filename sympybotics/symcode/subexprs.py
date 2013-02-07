
import sympy



class Subexprs(object):
       
  def __init__(self, mode='deep', ivars_name='tmp'):
    
    self.symbols = sympy.utilities.iterables.numbered_symbols(ivars_name, start=0, real=True)
    self.subexprs = list()
    if mode == 'deep':
      self._collect_func = self._collect_singleops
    elif mode == 'simple':
      self._collect_func = self._collect_nonatom
    
    
  def _collect_nonatom(self, expr):
    if expr.is_Atom:
      return expr
    else:
      ivar = next(self.symbols)
      self.subexprs.append((ivar,expr))
      return ivar
    
    
  def _add_if_new(self, expr):
    for (ivar,iexpr) in self.subexprs:
      if expr == iexpr:
        return ivar
    new_ivar = next(self.symbols)
    self.subexprs.append((new_ivar, expr))
    return new_ivar
  
  
  def _collect_singleops(self, expr):
    if expr.is_Atom:
        return expr
    else:
        new_args = []
        for arg in expr.args:
          new_args.append( self._collect_singleops(arg) )
        return self._add_if_new(type(expr)(*new_args))
    
    
  def collect(self, exprs):
    
    if isinstance(exprs, sympy.Basic): # if only one expression is passed
      exprs = [exprs]
      is_single_expr = True
    else:
      is_single_expr = False
    
    out_exprs = list()
    
    out_exprs.extend(map(self._collect_func, exprs))
          
    if is_single_expr:
      return out_exprs[0]
    elif isinstance(exprs, sympy.Matrix):
      return sympy.Matrix(exprs.rows, exprs.cols, out_exprs)
    else:
      return out_exprs
