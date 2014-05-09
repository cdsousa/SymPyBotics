
import subprocess


class VarMaps(object):
    pass

def gen_qepcad_varmaps( var_list ):

    forwmap = {}

    for v in var_list:

        v = str(v)
        if v in forwmap:
            raise 'Error: Two variables with same string representation'

        vm = v.replace('_','')

        if vm not in forwmap.values():
            forwmap[v] = vm
        else:
            raise 'Error: Two variables with same map'

    backmap = dict( zip( forwmap.values(), forwmap.keys()) )

    varmaps = VarMaps();

    varmaps.forward = forwmap
    varmaps.backward = backmap

    return varmaps



def sym_to_qepcad( symexpr, vardictforw=None ):
    if vardictforw is None:
        vardictforw = {}

    symexpr =  str(symexpr).replace('**','^').replace('*',' ')

    for var in vardictforw:
        symexpr = symexpr.replace(var,vardictforw[var])

    elements = symexpr.split()

    for i,e in enumerate(elements):
      if '/' in e:
        l,r = e.split('/')
        if not l.isnumeric():
          elements[i] = l + ' 1/' + r

    symexpr = ' '.join(elements)

    return symexpr



def qepcad_to_sym( qepcad_rel, vardictback=None ):
    if vardictforw is None:
        vardictforw = {}

    elements = qepcad_rel.split()

    out_elements = []

    last_isalnum = False

    for i,e in enumerate(elements):

        e0 = e.split('^')[0]
        e = e.replace('^','**')

        isalnum = e0.isalnum()

        if isalnum and last_isalnum:
            out_elements.append('*')

        last_isalnum = isalnum

        if isalnum and e0 in vardictback:
            e = e.replace( e0, str(vardictback[e0]) )

        out_elements.append( e )

    return ' '.join(out_elements)


def gen_qepcad_input( freevars, qvars, prenex, vardictforw=None ):
  if vardictforw is None:
      vardictforw = {}

  vars = '(' + ','.join( [ sym_to_qepcad(var) for var in freevars ] + [ sym_to_qepcad(var) for var in qvars ] ) + ')'

  for var in vardictforw:
    vars = vars.replace(var,vardictforw[var])

  freevarnum = str(len(freevars))

  prenex += '.' if prenex[-1] != '.' else ''


  qepcad_input = '[]\n'+vars+'\n'+freevarnum+'\n'+prenex+'\nfinish\n'

  return qepcad_input


def run_qepcad( qepcad_cmd, qepcad_input ):

    proc = subprocess.Popen(qepcad_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE )

    proc.stdin.write(qepcad_input.encode())

    out = proc.communicate()

    outstr = out[0].decode("utf-8")

    #print(outstr)

    outlines = outstr.splitlines()

    return outlines[ outlines.index('An equivalent quantifier-free formula:') + 2 ]
