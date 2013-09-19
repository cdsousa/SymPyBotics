import sympy
import numpy
from math import cos, sin

import sympybotics
from sympybotics._compatibility_ import exec_





def test_scara_dh_sym_geo_kin():

    pi = sympy.pi
    q = sympybotics.robotdef.q

    a1, a2, d3, d4 = sympy.symbols('a1, a2, d3, d4')

    scara = sympybotics.robotdef.RobotDef(
        'SCARA - Spong',
        [( 0, a1,  0, q),
         (pi, a2,  0, q),
         ( 0,  0,  q, 0),
         ( 0,  0, d4, q)],
        dh_convention='standard')

    scara_geo = sympybotics.geometry.Geometry(scara)
    scara_kin = sympybotics.kinematics.Kinematics(scara, scara_geo)

    cos, sin = sympy.cos, sympy.sin
    q1, q2, q3, q4 = sympy.flatten(scara.q)

    T_spong = sympy.Matrix([
        [(-sin(q1)*sin(q2) + cos(q1)*cos(q2))*cos(q4) + (sin(q1)*cos(q2) +
         sin(q2)*cos(q1))*sin(q4), -(-sin(q1)*sin(q2) +
         cos(q1)*cos(q2))*sin(q4) + (sin(q1)*cos(q2) +
         sin(q2)*cos(q1))*cos(q4), 0, a1*cos(q1) - a2*sin(q1)*sin(q2) +
         a2*cos(q1)*cos(q2)],
        [(sin(q1)*sin(q2) - cos(q1)*cos(q2))*sin(q4) + (sin(q1)*cos(q2) +
         sin(q2)*cos(q1))*cos(q4), (sin(q1)*sin(q2) -
         cos(q1)*cos(q2))*cos(q4) - (sin(q1)*cos(q2) +
         sin(q2)*cos(q1))*sin(q4), 0, a1*sin(q1) + a2*sin(q1)*cos(q2) +
         a2*sin(q2)*cos(q1)],
        [0, 0, -1, -d4 - q3],
        [0, 0, 0, 1]])

    J_spong = sympy.Matrix([[-a1*sin(q1) - a2*sin(q1)*cos(q2) -
                             a2*sin(q2)*cos(q1), -a2*sin(q1)*cos(q2) -
                             a2*sin(q2)*cos(q1), 0, 0],
                            [a1*cos(q1) - a2*sin(q1)*sin(q2) +
                             a2*cos(q1)*cos(q2), -a2*sin(q1)*sin(q2) +
                             a2*cos(q1)*cos(q2), 0, 0],
                            [0, 0, -1, 0],
                            [0, 0, 0, 0],
                            [0, 0, 0, 0],
                            [1, 1, 0, -1]])

    assert (scara_geo.T[-1] - T_spong).expand() == sympy.zeros(4)
    assert (scara_kin.J[-1] - J_spong).expand() == sympy.zeros((6, 4))


def test_puma_dh_num_geo_kin_dyn():
    pi = sympy.pi
    q = sympybotics.robotdef.q

    puma560_def = sympybotics.RobotDef(
        'Puma 560 Robot',
        [(   pi/2,        0,         0,   q),
         (      0,   0.4318,         0,   q),
         ('-pi/2', '0.0203', '0.15005', 'q'),  # test sympify
         (   pi/2,        0,    0.4318,   q),
         (  -pi/2,        0,         0,   q),
         (      0,        0,         0,   q)],
        dh_convention='standard')

    puma560_def.frictionmodel = None

    puma560 = sympybotics.RobotDynCode(puma560_def)

    q_test = [0.7504516826728697, 0.8395156106908136, 0.16851233582594916, 0.3849629637427072, 0.5252993946810777, 0.6701207256444748]
    dq_test = [0.24721855939629367, 0.9805915670454258, 0.9895299755642817, 0.7861135739668947, 0.273842245476577, 0.17182358900767503]
    ddq_test = [0.707405815485141, 0.25295715193420953, 0.9763909835998361, 0.8412822676113918, 0.4867768296473465, 0.11480270540937143]

    dynparm_test = [0.0, -0, -0, 0.34999999999999998, -0, 0.0, 0.0, 0.0, 0.0, 0, 1.03118515, 0.037980719999999996, 1.4401022999999999,
                    3.7274564059999999, -0.023751000000000001, 2.8425240559999998, -6.33012, 0.10439999999999999, 3.9585, 17.4,
                    0.090474288000000014, -0.0013739039999999998, 0.0068208000000000001, 0.111498032, 0.0047375999999999998,
                    0.015432319999999999, -0.09743999999999998, -0.06767999999999999, 0.336, 4.8, 0.0020960200000000001, -0, -0,
                    0.0012999999999999999, -0, 0.0020960200000000001, 0.0, 0.015579999999999998, 0.0, 0.82, 0.00029999999999999997, -0, -0,
                    0.00040000000000000002, -0, 0.00029999999999999997, 0.0, 0.0, 0.0, 0.34, 0.00024216, -0, -0, 0.00024216, -0,
                    4.0000000000000003e-05, 0.0, 0.0, 0.0028799999999999997, 0.09]

    dynparm_test2 = [0.47804562306292275, 0.5871876506259908, 0.9487349009746813, 0.35387185413632094, 0.28071604959871554, 0.368556182617345,
                     0.24355010647230801, 0.6753463418802456, 0.9728151452953864, 0.6620741264406734, 0.34669638996014096, 0.01593886435340608,
                     0.3521748260592804, 0.5384045845183812, 0.021600503502885116, 0.4654003203805651, 0.5202014161122065, 0.33744920539722967,
                     0.052363297799702835, 0.07051826001770234, 0.7389222546505236, 0.771434543548951, 0.3652539269897015, 0.2603059367721896,
                     0.4310648491411889, 0.7071252186366607, 0.16320122542325732, 0.44948421506462655, 0.48085540250421277, 0.08408482356412372,
                     0.923593157615906, 0.46852453511703684, 0.6670004526434297, 0.573634975268657, 0.12665814747855264, 0.3152822779549781,
                     0.21725524421221942, 0.5727451645276381, 0.7984025459787872, 0.7630923700903489, 0.27670044727992893, 0.919667661316697,
                     0.8828363383223908, 0.3618378476984413, 0.20690997158181246, 0.8321536435628591, 0.556153477173109, 0.3985225289975399,
                     0.7128303168401632, 0.433898343199971, 0.35116090526732047, 0.551636779637059, 0.03397585970570116, 0.28333059714010744,
                     0.4785548774382852, 0.16258401709779136, 0.8246056496848533, 0.23027686557606952, 0.26655111443159285, 0.6087446254045791]

    q_num_subs = dict(zip(puma560_def.q, q_test))

    T = puma560.geo.T[-1]
    T = T.subs(q_num_subs)
    T = numpy.matrix(T).astype(numpy.float64)

    J = puma560.kin.J[-1]
    J = J.subs(q_num_subs)
    J = numpy.matrix(J).astype(numpy.float64)

    M_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.M_code, 'M_puma560', puma560_def)
    c_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.c_code, 'c_puma560', puma560_def)
    g_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.g_code, 'g_puma560', puma560_def)
    tau_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.invdyn_code, 'tau_puma560', puma560_def)
    H_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.H_code, 'H_puma560', puma560_def)

    l = locals()
    exec_(M_func_def, globals(), l)
    exec_(c_func_def, globals(), l)
    exec_(g_func_def, globals(), l)
    exec_(tau_func_def, globals(), l)
    exec_(H_func_def, globals(), l)
    tau_puma560 = l['tau_puma560']
    g_puma560 = l['g_puma560']
    c_puma560 = l['c_puma560']
    M_puma560 = l['M_puma560']
    H_puma560 = l['H_puma560']

    tau = tau_puma560(dynparm_test, q_test, dq_test, ddq_test)
    tau = numpy.matrix(tau).T.astype(numpy.float64)

    g = g_puma560(dynparm_test, q_test)
    g = numpy.matrix(g).T.astype(numpy.float64)

    c = c_puma560(dynparm_test, q_test, dq_test)
    c = numpy.matrix(c).T.astype(numpy.float64)

    M = M_puma560(dynparm_test, q_test)
    M = numpy.matrix(M).reshape(puma560.dof, puma560.dof).astype(numpy.float64)

    H = H_puma560(q_test, dq_test, ddq_test)
    H = numpy.matrix(H).reshape(puma560.dof, puma560.dyn.n_dynparms
                                ).astype(numpy.float64)

    tau_t2 = tau_puma560(dynparm_test2, q_test, dq_test, ddq_test)
    tau_t2 = numpy.matrix(tau_t2).T.astype(numpy.float64)

    g_t2 = g_puma560(dynparm_test2, q_test)
    g_t2 = numpy.matrix(g_t2).T.astype(numpy.float64)

    c_t2 = c_puma560(dynparm_test2, q_test, dq_test)
    c_t2 = numpy.matrix(c_t2).T.astype(numpy.float64)

    M_t2 = M_puma560(dynparm_test2, q_test)
    M_t2 = numpy.matrix(M_t2).reshape(puma560.dof, puma560.dof
                                      ).astype(numpy.float64)

    T_pcorke = numpy.matrix([[-0.655113870655343, -0.474277925274361, -0.588120962109346, 0.0540498757911011],
                             [ 0.524340309498464,  0.275038163391149, -0.805866768463298, -0.154761559833806],
                             [ 0.543960528264717, -0.836310025215453, 0.0685017183295358,  0.568944734102513],
                             [               0.0,                0.0,                0.0,                1.0]])

    J_pcorke = numpy.matrix([[    0.154761559833806,   -0.416115317270121,   -0.181051500179202,                0.0,                0.0,                0.0],
                             [   0.0540498757911011,   -0.388002774727404,    -0.16881975145483,                0.0,                0.0,                0.0],
                             [-1.56125112837913e-17,  -0.0660115669805396,   -0.354377730397368,                0.0,                0.0,                0.0],
                             [ 2.77555756156289e-17,    0.681969181663071,    0.681969181663071, -0.618588326585383,  0.778592276991744, -0.588120962109346],
                             [ 1.04083408558608e-16,   -0.731380909828662,   -0.731380909828662, -0.576796808883883, -0.541217849732837, -0.805866768463298],
                             [                  1.0, 5.72458747072346e-17, 5.72458747072346e-17,  0.533529683779323,  0.317611878460765, 0.0685017183295358]])

    tau_pcorke = numpy.matrix([[  0.986185688341252],
                               [   16.1055550633721],
                               [  -6.54839494827661],
                               [0.00626452648190803],
                               [-0.0208972707132494],
                               [5.59805923645945e-5]])

    g_pcorke = numpy.matrix([[1.0295354054077e-15],
                             [   16.8135891283022],
                             [  -7.29033633567802],
                             [0.00449992182015992],
                             [ -0.026719893613902],
                             [                0.0]])

    c_pcorke = numpy.matrix([[   -0.240237883763956],
                             [     -1.0843848843681],
                             [    0.365583586254475],
                             [-0.000123404998742911],
                             [  0.00359951257553373],
                             [   1.1075729558421e-5]])

    M_pcorke = numpy.matrix([[    2.05591694561963,    -0.607149114401804,   -0.0772414386669633,   0.00110264316988551, 0.000264135330891356,  2.74006873318143e-6],
                             [  -0.607149114401804,      2.00048563124708,     0.306870038373922, -0.000227846077158843, 0.000780985328855223,  7.53260796283046e-6],
                             [ -0.0772414386669632,     0.306870038373922,     0.361368447500762, -0.000267122764221431,  0.00156301630001611,  7.53260796283046e-6],
                             [ 0.00110264316988551, -0.000227846077158844, -0.000267122764221432,   0.00169083802882858, 1.45380752747483e-20,   3.4606953693396e-5],
                             [0.000264135330891356,  0.000780985328855223,   0.00156301630001611,  2.33135442795155e-20,           0.00064216, 2.44929359829471e-21],
                             [ 2.74006873318143e-6,   7.53260796283046e-6,   7.53260796283046e-6,    3.4606953693396e-5, 2.44929359829471e-21,               4.0e-5]])

    tau_assrt2 = numpy.matrix([[ 4.55384377546122],
                               [-18.3482765770679],
                               [-18.5150406032816],
                               [-3.48354572715293],
                               [-4.22434630683546],
                               [-8.50580534103007]])

    g_assrt2 = numpy.matrix([[-4.44089209850063e-16],
                             [    -16.5628257742538],
                             [    -23.0524142097215],
                             [     -4.2536826692411],
                             [    -8.37451830056579],
                             [     -7.9940463468124]])

    c_assrt2 = numpy.matrix([[   3.79031173486171],
                             [  -6.37866680030809],
                             [  0.668430906458299],
                             [ 0.0446113764650349],
                             [   2.68800462309735],
                             [-0.0871119983741941]])

    M_assrt2 = numpy.matrix([[  2.80338554010218, -0.884144436072495,  -1.36196693847736,   0.914512998462046, -0.799840112619873,  -0.402048288616148],
                             [-0.884144436072495,   5.09045072475534,   3.75033319331688,  -0.822003775709123,   2.00576876282064,  -0.136034063369629],
                             [ -1.36196693847736,   3.75033319331688,   3.69036945486662,  -0.784876338891255,   1.91719457979014,  0.0657267986145439],
                             [ 0.914512998462046, -0.822003775709123, -0.784876338891255,    1.86573777350968,  -1.06272697183618, 0.00496876753726561],
                             [-0.799840112619873,   2.00576876282064,   1.91719457979014,   -1.06272697183618,    1.2083737144556,  -0.396167549434339],
                             [-0.402048288616148, -0.136034063369629, 0.0657267986145439, 0.00496876753726561, -0.396167549434339,   0.162584017097791]])

    assert_precision = 10

    assert not numpy.any(numpy.round(T_pcorke - T, assert_precision))
    assert not numpy.any(numpy.round(J_pcorke - J, assert_precision))
    assert not numpy.any(numpy.round(tau_pcorke - tau, assert_precision))
    assert not numpy.any(numpy.round(g_pcorke - g, assert_precision))
    assert not numpy.any(numpy.round(c_pcorke - c, assert_precision))
    assert not numpy.any(numpy.round(M_pcorke - M, assert_precision))
    assert not numpy.any(numpy.round(
        tau_pcorke - H * numpy.matrix(dynparm_test).T, assert_precision))
    assert not numpy.any(numpy.round(tau_assrt2 - tau_t2, assert_precision))
    assert not numpy.any(numpy.round(g_assrt2 - g_t2, assert_precision))
    assert not numpy.any(numpy.round(c_assrt2 - c_t2, assert_precision))
    assert not numpy.any(numpy.round(M_assrt2 - M_t2, assert_precision))
    assert not numpy.any(numpy.round(
        tau_assrt2 - H * numpy.matrix(dynparm_test2).T, assert_precision))


def test_puma_mdh_num_geo_kin_dyn():
    pi = sympy.pi
    q = sympybotics.robotdef.q

    puma560_def = sympybotics.RobotDef(
        'Puma Robot 560 mdh',
        [(    0.,      0.,      0., q),
         (-pi/2.,      0.,  0.2435, q),
         (    0.,  0.4318, -0.0934, q),
         ( pi/2., -0.0203,  0.4331, q),
         (-pi/2.,      0.,      0., q),
         ( pi/2.,      0.,      0., q)],
        dh_convention='modified')

    puma560_def.frictionmodel = None

    puma560 = sympybotics.RobotDynCode(puma560_def)

    q_test = [0.7504516826728697, 0.8395156106908136, 0.16851233582594916, 0.3849629637427072, 0.5252993946810777, 0.6701207256444748]
    dq_test = [0.24721855939629367, 0.9805915670454258, 0.9895299755642817, 0.7861135739668947, 0.273842245476577, 0.17182358900767503]
    ddq_test = [0.707405815485141, 0.25295715193420953, 0.9763909835998361, 0.8412822676113918, 0.4867768296473465, 0.11480270540937143]

    dynparm_test = [0.0, 0.0, 0.0, 0.0, 0.0, 0.34999999999999998, 0, 0, 0, 0, 0.1350808, -0.0070992, 0.018931199999999999, 0.60891200000000001,
                    0.0016703999999999998, 0.62008400000000008, 1.1832, 0.10439999999999999, -0.2784, 17.4, 0.090460800000000008, 0.0, 0.0,
                    0.013440800000000001, 0.0047039999999999998, 0.089520000000000002, 0.0, -0.336, 0.0672, 4.8, 0.0020960200000000001, 0.0,
                    0.0, 0.0020960200000000001, 0.0, 0.0012999999999999999, 0.0, 0.0, -0.015579999999999998, 0.82, 0.00029999999999999997, 0.0,
                    0.0, 0.00029999999999999997, 0.0, 0.00040000000000000002, 0.0, 0.0, 0.0, 0.34, 0.00024216, 0.0, 0.0, 0.00024216, 0.0,
                    4.0000000000000003e-05, 0.0, 0, 0.0028799999999999997, 0.09]

    dynparm_test2 = [0.47804562306292275, 0.5871876506259908, 0.9487349009746813, 0.35387185413632094, 0.28071604959871554, 0.368556182617345,
                     0.24355010647230801, 0.6753463418802456, 0.9728151452953864, 0.6620741264406734, 0.34669638996014096, 0.01593886435340608,
                     0.3521748260592804, 0.5384045845183812, 0.021600503502885116, 0.4654003203805651, 0.5202014161122065, 0.33744920539722967,
                     0.052363297799702835, 0.07051826001770234, 0.7389222546505236, 0.771434543548951, 0.3652539269897015, 0.2603059367721896,
                     0.4310648491411889, 0.7071252186366607, 0.16320122542325732, 0.44948421506462655, 0.48085540250421277, 0.08408482356412372,
                     0.923593157615906, 0.46852453511703684, 0.6670004526434297, 0.573634975268657, 0.12665814747855264, 0.3152822779549781,
                     0.21725524421221942, 0.5727451645276381, 0.7984025459787872, 0.7630923700903489, 0.27670044727992893, 0.919667661316697,
                     0.8828363383223908, 0.3618378476984413, 0.20690997158181246, 0.8321536435628591, 0.556153477173109, 0.3985225289975399,
                     0.7128303168401632, 0.433898343199971, 0.35116090526732047, 0.551636779637059, 0.03397585970570116, 0.28333059714010744,
                     0.4785548774382852, 0.16258401709779136, 0.8246056496848533, 0.23027686557606952, 0.26655111443159285, 0.6087446254045791]

    q_num_subs = dict(zip(puma560_def.q, q_test))

    T = puma560.geo.T[-1]
    T = T.subs(q_num_subs)
    T = numpy.matrix(T).astype(numpy.float64)

    J = puma560.kin.J[-1]
    J = J.subs(q_num_subs)
    J = numpy.matrix(J).astype(numpy.float64)

    M_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.M_code, 'M_puma560', puma560_def)
    c_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.c_code, 'c_puma560', puma560_def)
    g_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.g_code, 'g_puma560', puma560_def)
    tau_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.invdyn_code, 'tau_puma560', puma560_def)
    H_func_def = sympybotics.robotcodegen.dyn_code_to_func(
        'python', puma560.H_code, 'H_puma560', puma560_def)

    l = locals()
    exec_(M_func_def, globals(), l)
    exec_(c_func_def, globals(), l)
    exec_(g_func_def, globals(), l)
    exec_(tau_func_def, globals(), l)
    exec_(H_func_def, globals(), l)
    tau_puma560 = l['tau_puma560']
    g_puma560 = l['g_puma560']
    c_puma560 = l['c_puma560']
    M_puma560 = l['M_puma560']
    H_puma560 = l['H_puma560']

    tau = tau_puma560(dynparm_test, q_test, dq_test, ddq_test)
    tau = numpy.matrix(tau).T.astype(numpy.float64)

    g = g_puma560(dynparm_test, q_test)
    g = numpy.matrix(g).T.astype(numpy.float64)

    c = c_puma560(dynparm_test, q_test, dq_test)
    c = numpy.matrix(c).T.astype(numpy.float64)

    M = M_puma560(dynparm_test, q_test)
    M = numpy.matrix(M).reshape(puma560.dof, puma560.dof).astype(numpy.float64)

    H = H_puma560(q_test, dq_test, ddq_test)
    H = numpy.matrix(H).reshape(puma560.dof, puma560.dyn.n_dynparms
                                ).astype(numpy.float64)

    tau_t2 = tau_puma560(dynparm_test2, q_test, dq_test, ddq_test)
    tau_t2 = numpy.matrix(tau_t2).T.astype(numpy.float64)

    g_t2 = g_puma560(dynparm_test2, q_test)
    g_t2 = numpy.matrix(g_t2).T.astype(numpy.float64)

    c_t2 = c_puma560(dynparm_test2, q_test, dq_test)
    c_t2 = numpy.matrix(c_t2).T.astype(numpy.float64)

    M_t2 = M_puma560(dynparm_test2, q_test)
    M_t2 = numpy.matrix(M_t2).reshape(puma560.dof, puma560.dof
                                      ).astype(numpy.float64)

    T_pcorke = numpy.matrix([[-0.655113870655343, -0.474277925274361,  0.588120962109346,   0.368531204501578],
                             [ 0.524340309498464,  0.275038163391149,  0.805866768463298,   0.548861637700235],
                             [-0.543960528264717,  0.836310025215453, 0.0685017183295358, -0.0731561881633544],
                             [               0.0,                0.0,                0.0,                 1.0]])

    J_pcorke = numpy.matrix([[   -0.548861637700235,   -0.053505039458511,    0.181558777632408,               0.0,                0.0,                0.0],
                             [    0.368531204501578,  -0.0498902657753524,    0.169292757497222,               0.0,                0.0,                0.0],
                             [ 2.66713734431434e-17,   -0.643843409557299,    -0.35547724614047,               0.0,                0.0,                0.0],
                             [-2.77555756156289e-17,   -0.681969181663072,   -0.681969181663072, 0.618588326585383, -0.778592276991744,  0.588120962109346],
                             [-1.04083408558608e-16,    0.731380909828662,    0.731380909828662, 0.576796808883883,  0.541217849732837,  0.805866768463298],
                             [                  1.0, 5.72458747072346e-17, 5.72458747072346e-17, 0.533529683779323,  0.317611878460765, 0.0685017183295358]])

    tau_pcorke = numpy.matrix([[   2.71610302653938],
                               [  -28.7944111980059],
                               [  -7.07013200484716],
                               [0.00682888994276653],
                               [-0.0220889100229662],
                               [5.59805923645946e-5]])

    g_pcorke = numpy.matrix([[2.88527195030991e-15],
                             [   -31.1492065095002],
                             [   -7.04528104551999],
                             [ 0.00449992182016002],
                             [  -0.026719893613902],
                             [                 0.0]])

    c_pcorke = numpy.matrix([[   0.634177626441184],
                             [    1.12398760139525],
                             [  -0.373925639216401],
                             [-8.54577359118394e-5],
                             [ 0.00162606594660281],
                             [  1.1075729558421e-5]])

    M_pcorke = numpy.matrix([[   2.88092539107775,     0.454967124647697,   -0.0748835030743064,   0.00187610810095465,  0.00080480941146298,  2.74006873318143e-6],
                             [  0.454967124647696,      2.16866309038072,     0.368189596401607, -0.000307104503389342,  0.00234804941178354,  7.53260796283046e-6],
                             [-0.0748835030743065,     0.368189596401607,     0.315830104422493, -0.000267827816326753,  0.00156601844062265,  7.53260796283046e-6],
                             [0.00187610810095465, -0.000307104503389342, -0.000267827816326753,   0.00169083802882858, 1.98357371583556e-20,   3.4606953693396e-5],
                             [0.00080480941146298,   0.00234804941178354,   0.00156601844062265,  2.86112061631228e-20,           0.00064216, 2.44929359829471e-21],
                             [2.74006873318143e-6,   7.53260796283046e-6,   7.53260796283046e-6,    3.4606953693396e-5, 2.44929359829471e-21,               4.0e-5]])

    tau_assrt2 = numpy.matrix([[  5.46017281703168],
                               [ -8.53813836226581],
                               [-0.490552199963579],
                               [  16.4258371315422],
                               [  6.18803453956645],
                               [  7.99017179545593]])

    g_assrt2 = numpy.matrix([[-1.77635683940025e-15],
                             [    -12.2719368424473],
                             [    -5.98183554268347],
                             [     18.5386690234033],
                             [    0.512467781908107],
                             [      7.9940463468124]])

    c_assrt2 = numpy.matrix([[  3.69185531689127],
                             [ -1.34087787458136],
                             [  1.44951875057327],
                             [ -2.77730540324791],
                             [  3.31210163468715],
                             [-0.788977149763209]])

    M_assrt2 = numpy.matrix([[  3.65840168922149, -0.555774017452474,    -1.39220275090878,  0.985927680308992, -0.439435427090482,    0.563816026114764],
                             [-0.555774017452474,   5.33926241554461,     4.32747701210166,  -1.24139932228482,   1.97231403813762,   -0.207221260322988],
                             [ -1.39220275090878,   4.32747701210166,     4.13345124166948,  -1.07958817736071,   1.65363538346271, -0.00546039833881495],
                             [ 0.985927680308992,  -1.24139932228482,    -1.07958817736071,   1.75689881264024, -0.358749490679038,    0.276358110012213],
                             [-0.439435427090482,   1.97231403813762,     1.65363538346271, -0.358749490679038,   1.67868951032001,    0.396167549434339],
                             [ 0.563816026114764, -0.207221260322988, -0.00546039833881495,  0.276358110012213,  0.396167549434339,    0.162584017097791]])

    assert_precision = 10

    assert not numpy.any(numpy.round(T_pcorke - T, assert_precision))
    assert not numpy.any(numpy.round(J_pcorke - J, assert_precision))
    assert not numpy.any(numpy.round(tau_pcorke - tau, assert_precision))
    assert not numpy.any(numpy.round(g_pcorke - g, assert_precision))
    assert not numpy.any(numpy.round(c_pcorke - c, assert_precision))
    assert not numpy.any(numpy.round(M_pcorke - M, assert_precision))
    assert not numpy.any(numpy.round(
        tau_pcorke - H * numpy.matrix(dynparm_test).T, assert_precision))
    assert not numpy.any(numpy.round(tau_assrt2 - tau_t2, assert_precision))
    assert not numpy.any(numpy.round(g_assrt2 - g_t2, assert_precision))
    assert not numpy.any(numpy.round(c_assrt2 - c_t2, assert_precision))
    assert not numpy.any(numpy.round(M_assrt2 - M_t2, assert_precision))
    assert not numpy.any(numpy.round(
        tau_assrt2 - H * numpy.matrix(dynparm_test2).T, assert_precision))
