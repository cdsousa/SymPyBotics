""" Extra dynamic terms modelation """

from sympy import zeros, sign
from ..utils import identity


_frictionterms = set(['Coulomb', 'viscous', 'offset'])


def frictionforce(rbtdef, ifunc=None):
    '''Generate friction forces (Coulomb/viscouse model plus offset).'''

    if not ifunc:
        ifunc = identity

    fric = zeros(rbtdef.dof, 1)

    if rbtdef.frictionmodel is None or len(rbtdef.frictionmodel) == 0:
        pass
    else:
        askedterms = set(rbtdef.frictionmodel)
        if askedterms.issubset(_frictionterms):
            if 'viscous' in askedterms:
                for i in range(rbtdef.dof):
                    fric[i] += rbtdef.fv[i] * rbtdef.dq[i]
            if 'Coulomb' in askedterms:
                for i in range(rbtdef.dof):
                    fric[i] += rbtdef.fc[i] * sign(rbtdef.dq[i])
            if 'offset' in askedterms:
                for i in range(rbtdef.dof):
                    fric[i] += rbtdef.fo[i]
            fric[i] = ifunc(fric[i])
        else:
            raise Exception(
                'Friction model terms \'%s\' not understanded. Use None or a'
                ' combination of %s.' %
                (str(askedterms - _frictionterms), _frictionterms))

    return fric


def driveinertiaterm(rbtdef, ifunc=None):
    '''Generate drive inertia term (Siplified, neglets gyroscopic effects).'''

    if not ifunc:
        ifunc = identity

    driveinertia = zeros(rbtdef.dof, 1)

    if rbtdef.driveinertiamodel is None:
        pass
    elif rbtdef.driveinertiamodel == 'simplified':
        for i in range(rbtdef.dof):
            driveinertia[i] = rbtdef.Ia[i] * rbtdef.ddq[i]
            driveinertia[i] = ifunc(driveinertia[i])
    else:
        raise Exception('Drive inertia model \'%s\' not understanded. Use'
                        ' None or \'simplified\'.' % rbtdef.driveinertiamodel)

    return driveinertia
