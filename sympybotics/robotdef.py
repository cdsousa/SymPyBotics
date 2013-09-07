
import sympy


def _new_sym(name):
    return sympy.symbols(name, real=True)


def _elements_to_tensor(elems):
    return sympy.Matrix([[elems[0], elems[1], elems[2]],
                         [elems[1], elems[3], elems[4]],
                         [elems[2], elems[4], elems[5]]])


def _elementslist_to_tensorlist(elementslist):

    return [_elements_to_tensor(elems) for elems in elementslist]


def _sym_skew(v):
    return sympy.Matrix([[0, -v[2],  v[1]],
                         [v[2],     0, -v[0]],
                         [-v[1],  v[0],     0]])

_joint_symb = _new_sym('q')


def _joint_i_symb(i):
    return _new_sym('q' + str(i))

q = _joint_symb

_dh_alpha, _dh_a, _dh_d, _dh_theta = _new_sym('alpha,a,d,theta')
default_dh_symbols = (_dh_alpha, _dh_a, _dh_d, _dh_theta)

_cos = sympy.cos
_sin = sympy.sin
_standard_dh_transfmat = sympy.Matrix([
    [_cos(_dh_theta), -_sin(_dh_theta) * _cos(_dh_alpha),
     _sin(_dh_theta) * _sin(_dh_alpha), _dh_a * _cos(_dh_theta)],
    [_sin(_dh_theta),  _cos(_dh_theta) * _cos(_dh_alpha), -
     _cos(_dh_theta) * _sin(_dh_alpha), _dh_a * _sin(_dh_theta)],
    [0, _sin(_dh_alpha), _cos(_dh_alpha), _dh_d],
    [0, 0, 0, 1]])
_modified_dh_transfmat = sympy.Matrix([
    [_cos(_dh_theta), -_sin(_dh_theta), 0, _dh_a],
    [_sin(_dh_theta) * _cos(_dh_alpha), _cos(_dh_theta) *
     _cos(_dh_alpha), -_sin(_dh_alpha), -_sin(_dh_alpha) * _dh_d],
    [_sin(_dh_theta) * _sin(_dh_alpha), _cos(_dh_theta) *
     _sin(_dh_alpha), _cos(_dh_alpha), _cos(_dh_alpha) * _dh_d],
    [0, 0, 0, 1]])


class RobotDef(object):

    """Class that generates and holds robot definitions and symbols."""

    def __init__(self, name, dh_parms, dh_convention, shortname=None):
        """
        Create RobotDef instance with data structures for robot geometry and
        symbols of robot dynamics.
        """

        self.dof = len(dh_parms)

        self.name = str(name)
        if shortname is not None:
            self.shortname = str(shortname).replace(' ', '_').replace('.', '_')
        else:
            self.shortname = ''.join(c for c in name.replace(
                ' ', '_').replace('.', '_') if c.isalnum() or c == '_')

        dh_convention = dh_convention.lower()
        if dh_convention in ['standard', 'std', 'dh', 'sdh']:
            self._dh_convention = 'standard'
            self._dh_transfmat = _standard_dh_transfmat
        elif dh_convention in ['modified', 'mod', 'mdh']:
            self._dh_convention = 'modified'
            self._dh_transfmat = _modified_dh_transfmat
        else:
            raise ValueError(
                "DH convention %s not known/implemented (use 'standard' or"
                "'modified')" % dh_convention)

        self._dh_symbols = default_dh_symbols

        self._dyn_parms_order = 'Khalil'

        self.frictionmodel = None  # can be None or 'simple'

        # g_a = sympy.symbols('g_a',real=True)
        g_a = 9.81
        self.gravity = sympy.Matrix([[0.0], [0.0], [-g_a]])

        self._gen_symbols()

        self._set_dh_parms(dh_parms)

    @property
    def dh_convention(self):
        return self._dh_convention

    @property
    def dyn_parms_order(self):
        return self._dyn_parms_order

    @property
    def dh_parms(self):
        return self._dh_parms

    def __str__(self):
        return 'RobotDef instance: ' + self.name

    @property
    def L(self):
        return _elementslist_to_tensorlist(self.Le)

    @property
    def I(self):
        return _elementslist_to_tensorlist(self.Ie)

    def _gen_symbols(self):
        """Generate robot dynamic symbols and populates RobotDef instance with
        them. (internal function)"""

        dof = self.dof

        # q = self.q = [ _new_sym('q_'+str(i+1)) for i in range(self.dof) ]
        # dq = self.dq = [ _new_sym('dq_'+str(i+1)) for i in range(self.dof) ]
        # ddq = self.ddq = [ _new_sym('ddq_'+str(i+1)) for i in range(self.dof)
        # ]
        self.q = sympy.Matrix([[_new_sym('q' + str(i + 1))]
                              for i in range(self.dof)])
        self.dq = sympy.Matrix([[_new_sym(r'\dot{q}_' + str(i + 1))]
                               for i in range(self.dof)])
        self.ddq = sympy.Matrix(
            [[_new_sym(r'\ddot{q}_' + str(i + 1))] for i in range(self.dof)])

        self.non_latex_symbols = {}
        for i in range(self.dof):
            self.non_latex_symbols[self.q[i]] = r'q' + str(i + 1)
            self.non_latex_symbols[
                self.dq[i]] = r'dq' + str(i + 1)
            self.non_latex_symbols[
                self.ddq[i]] = r'ddq' + str(i + 1)

        m = self.m = list(range(self.dof))
        l = self.l = list(range(self.dof))
        Le = self.Le = list(range(self.dof))

        r = self.r = list(range(self.dof))
        Ie = self.Ie = list(range(self.dof))

        fv = self.fv = list(range(self.dof))
        fc = self.fc = list(range(self.dof))

        for i in range(dof):

            m[i] = _new_sym('m_' + str(i + 1))
            l[i] = sympy.Matrix([_new_sym('l_' + str(i + 1) + dim)
                                for dim in ['x', 'y', 'z']])
            Le[i] = [_new_sym('L_' + str(i + 1) + elem)
                     for elem in ['xx', 'xy', 'xz', 'yy', 'yz', 'zz']]

            r[i] = sympy.Matrix([_new_sym('r_' + str(i + 1) + dim)
                                for dim in ['x', 'y', 'z']])
            Ie[i] = [_new_sym('I_' + str(i + 1) + elem)
                     for elem in ['xx', 'xy', 'xz', 'yy', 'yz', 'zz']]

            fv[i] = _new_sym('fv_' + str(i + 1))
            fc[i] = _new_sym('fc_' + str(i + 1))

        I = self.I
        L = self.L

        I_funcof_L = self.I_funcof_L = list(range(self.dof))
        L_funcof_I = self.L_funcof_I = list(range(self.dof))

        dict_I2Lexp = self.dict_I2Lexp = dict()
        dict_L2Iexp = self.dict_L2Iexp = dict()
        dict_l2mr = self.dict_l2mr = dict()
        dict_r2lm = self.dict_r2lm = dict()

        for i in range(dof):

            L_funcof_I[i] = I[i] + m[i] * _sym_skew(r[i]).T * _sym_skew(r[i])
            I_funcof_L[i] = L[i] - m[i] * _sym_skew(r[i]).T * _sym_skew(r[i])

            for elem, exprss in enumerate(I_funcof_L[i]):
                dict_I2Lexp[I[i][elem]] = exprss

            for elem, exprss in enumerate(L_funcof_I[i]):
                dict_L2Iexp[L[i][elem]] = exprss

            for elem in range(3):
                dict_l2mr[l[i][elem]] = m[i] * r[i][elem]
                dict_r2lm[r[i][elem]] = l[i][elem] / m[i]

        return self

    def _set_dh_parms(self, dh_parms_list):
        """
        Define the RobotDef geometry using Denavit-Hartenberg notation.
        """

        if len(dh_parms_list) != self.dof:
            raise Exception('RobotDef.set_geometry(): provided number of links'
                            'differ from robot dof (%d vs %d).' %
                            (len(dh_parms_list), self.dof))

        self._dh_parms = []

        self._links_sigma = [0] * self.dof
        theta_index = self._dh_symbols.index(sympy.Symbol('theta', real=True))
        d_index = self._dh_symbols.index(sympy.Symbol('d', real=True))

        for i in range(self.dof):

            if len(dh_parms_list[i]) != 4:
                raise Exception(
                    'RobotDef.set_dh_parms: wrong number of Denavit-Hartenberg'
                    'parameters (must be 4 per link).')

            link_dh_parms = []
            for p in dh_parms_list[i]:

                p = sympy.sympify(p)

                for s in p.free_symbols:
                    if str(s) == str(_joint_symb):
                        p = p.subs(s, _joint_i_symb(i + 1))

                link_dh_parms.append(p)

            self._dh_parms.append(tuple(link_dh_parms))

            try:
                if self._dh_parms[i][theta_index].has(self.q[i]):
                    self._links_sigma[i] = 0
                    # print 'joint',i+1,'is revolute'
            except:
                pass
            try:
                if self._dh_parms[i][d_index].has(self.q[i]):
                    self._links_sigma[i] = 1
                    # print 'joint',il+1,'is prismatic'
            except:
                pass

        return self

    def dynparms(self, parm_order=None):
        """Return list of RobotDef symbolic dynamic parameters."""

        if not parm_order:
            parm_order = self._dyn_parms_order
        parm_order = parm_order.lower()

        parms = []
        for i in range(0, self.dof):

            if parm_order == 'khalil' or parm_order == 'tensor first':
                # Lxx Lxy Lxz Lyy Lyz Lzz lx ly lz m
                parms += self.Le[i]
                parms += sympy.flatten(self.l[i])
                parms += [self.m[i]]

            elif parm_order == 'siciliano' or parm_order == 'mass first':
                # m lx ly lz Lxx Lxy Lxz Lyy Lyz Lzz
                parms += [self.m[i]]
                parms += sympy.flatten(self.l[i])
                parms += self.Le[i]

            else:
                raise Exception(
                    'RobotDef.Parms(): dynamic parameters order \''
                    + parm_order + '\' not know.')

            if self.frictionmodel == 'simple':
                parms += [self.fv[i], self.fc[i]]

        return parms

    def _set_new_dynparms(self, dynparms):
        """Define new symbols to dynamic parameters."""
        symparms = self.dynparms()
        parmsidx = dict(zip(symparms, range(len(symparms))))
        for i in range(self.dof):
            for x, Le_xx in enumerate(self.Le[i]):
                self.Le[i][x] = dynparms[parmsidx[Le_xx]]
            for x, l_x in enumerate(self.l[i]):
                self.l[i][x] = dynparms[parmsidx[l_x]]
            self.m[i] = dynparms[parmsidx[self.m[i]]]
            if self.frictionmodel == 'simple':
                self.fv[i] = dynparms[parmsidx[self.fv[i]]]
                self.fc[i] = dynparms[parmsidx[self.fc[i]]]
