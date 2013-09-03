
class NotAvailableError(Exception):
    def __init__(self, function_name):
        msg = 'Function %s not available since cvxopt package '\
              'was not found' % function_name
        Exception.__init__(self, msg)

try:
    import cvxopt
except ImportError:
    cvxopt = None
else:
    import cvxopt
    import cvxopt.solvers


import numpy
import sympy


def mrepl(m,repl):
  return m.applyfunc(lambda x: x.xreplace(repl))

def skew(v):
  return sympy.Matrix( [ [     0, -v[2],  v[1] ],
                         [  v[2],     0, -v[0] ],
                         [ -v[1],  v[0],     0 ] ] )

def regr_matrices( dof, parm_num, q, dq, ddq, tau, regr_func):

  sn = q.shape[0]

  H_S = numpy.matrix( numpy.zeros( ( dof*sn, parm_num ) ) )
  tau_S = numpy.matrix( numpy.zeros( dof*sn ) ).T

  for i in range(sn):
    H_S[ i*dof : i*dof+dof , : ] = numpy.array( regr_func( q[i], dq[i], ddq[i] ) ).reshape(dof,parm_num)

  for i in range(sn):
      tau_S[ i*dof : i*dof+dof ] = numpy.mat( tau[i] ).T

  return H_S,tau_S




def get_diag_blocks(M):
  b = []
  n = M.shape[0]
  c = 0
  while c < n-1:
      for l in range(n-1,c,-1):
          if M[l,c] != 0 or M[c,l] != 0:
              break
          elif l == c+1:
              b.append(l)
      c = l
  return b + [n]





def prepare_sdp( var_symbs, LMI_matrix, split_diag_blocks=True):
    v = var_symbs
    vn = len(v)

    if isinstance( LMI_matrix, list ) or isinstance( LMI_matrix, tuple ):

      blocks_LMI_matrix = LMI_matrix

    else:

      if split_diag_blocks:

          blocks = get_diag_blocks(LMI_matrix)
          blocks_LMI_matrix = []
          for a,b in zip( [0]+blocks[:-1] , blocks ):
              blocks_LMI_matrix.append( LMI_matrix[ a:b, a:b ] )
          print('Split into %d diagonal blocks.'%(len(blocks)))

      else:

          blocks_LMI_matrix = [ LMI_matrix ]


    blocks_Fi = []

    zero_subs = dict(zip(v,[0]*vn))

    for LMI_matrix in blocks_LMI_matrix:

        Fi = [0]*(1+vn)

        Fsym = mrepl( LMI_matrix, zero_subs )
        Fi[0] = numpy.matrix(Fsym).astype(float)

        for i,s in enumerate(v):
            Fsym = LMI_matrix.applyfunc(lambda x: 0 if not x.coeff(s) else x.coeff(s))
            Fi[i+1] = numpy.matrix(Fsym).astype(float)

        blocks_Fi.append( Fi )

    return blocks_Fi


def sdp( c, Ai_blocks, primalstar=None, dualstart=None, solver='dsdp', verbose=0, interpret=False, maxiters=1000, tolerance=10e-7 ):
  if cvxopt is None:
    raise NotAvailableError(sdp.__name__)


  c = cvxopt.matrix( c )

  h = []
  G = []

  for A in Ai_blocks:
      h.append( cvxopt.matrix( A[0].astype(float).tolist() ) )
      G.append( cvxopt.matrix( [ (-A[i]).flatten().astype(float).tolist()[0] for i in range(1,len(A)) ] ) )

  if solver != 'dsdp' and solver != 'conelp':
      raise Exception('error: unknown solver (available: dsdp or conelp)')

  if solver == 'conelp': solver == ''

  cvxopt.solvers.options['show_progress'] = (1 if verbose > 0 else 0) #True/False (default: True)
  cvxopt.solvers.options['maxiters'] = maxiters #positive integer (default: 100)
  # cvxopt.solvers.options['refinement'] #positive integer (default: 1)
  cvxopt.solvers.options['abstol'] = tolerance #scalar (default: 1e-7)
  cvxopt.solvers.options['reltol'] = tolerance #scalar (default: 1e-6)
  cvxopt.solvers.options['feastol'] = tolerance #scalar (default: 1e-7).

  cvxopt.solvers.options['DSDP_Monitor'] = (1 if verbose > 0 else 0) #True/False (default: 0)
  cvxopt.solvers.options['DSDP_MaxIts'] = maxiters #positive integer
  cvxopt.solvers.options['DSDP_GapTolerance'] = tolerance #scalar (default: 1e-5).

  sol = cvxopt.solvers.sdp(c, Gs=G, hs=h, solver=solver)

  if not interpret:
      if verbose >= 0: print(sol['status'])
      return sol

  if verbose > 0:
      print('------------------------------\nsol[\'status\'] : '+sol['status']+'\n\n____________________________________')

  if sol['status'] == 'optimal':
      if verbose >= 0: print('optimal')
      return numpy.matrix(sol['x'])

  else:

      sp = 'residual as primal infeasibility certificate'
      sd = 'residual as dual infeasibility certificate'

      if solver == 'dsdp' and sol['status'] == 'unknown':

          if sol[sp] is not None:
              s = 'primal infeasible'
          elif sol[sp] is not None:
              s = 'dual infeasible'
          else:
              s = 'unknown'

          if verbose >= 0: print(s)
          if verbose >= 0: print(sp + ': ' + str(sol[sp]) + '\n' + sd + ': ' + str(sol[sd]))
          return


      else:

          if verbose >= 0: print(sol['status'])
          if verbose >= 0: print(sp + ': ' + str(sol[sp]) + '\n' + sd + ': ' + str(sol[sd]))
          return


