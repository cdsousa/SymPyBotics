
import sympy



class Subexprs(object):
       
  def __init__(self, mode='unique_ops', ivars_name='tmp'):
    
    self.symbols = sympy.utilities.iterables.numbered_symbols(ivars_name, start=0, real=True)
    self.subexprs = list()
    self.exprs_dict = dict()
    
    if mode == 'unique_ops':
      self._collect_func = self._collect_uniqueops
    elif mode == 'whole_exprs':
      self._collect_func = self._collect_exprs
    else:
      raise Exception("No '%s' sub-expressions collection mode known."%mode)
    
    
  def _collect_exprs(self, expr):
    if expr.is_Atom:
      return expr
    else:
      ivar = next(self.symbols)
      self.subexprs.append((ivar,expr))
      return ivar
    
    
  def _collect_op(self, expr):
    if expr in self.exprs_dict:
      return self.exprs_dict[expr]
    else:
      new_ivar = next(self.symbols)
      self.subexprs.append((new_ivar, expr))
      self.exprs_dict[expr] = new_ivar
      return new_ivar
  
  
  def _collect_uniqueops(self, expr):
    if expr.is_Atom:
        return expr
    else:
        new_args = []
        for arg in expr.args:
          new_args.append( self._collect_uniqueops(arg) )
        return self._collect_op(type(expr)(*new_args))
    
    
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
