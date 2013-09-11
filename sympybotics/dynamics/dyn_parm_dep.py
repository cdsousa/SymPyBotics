import numpy


def find_dyn_parm_deps(dof, parm_num, regressor_func):
    '''
    Find dynamic parameter dependencies (i.e., regressor column dependencies).
    '''

    samples = 10000
    round = 10

    pi = numpy.pi

    Z = numpy.zeros((dof * samples, parm_num))

    for i in range(samples):
        q = [float(numpy.random.random() * 2.0 * pi - pi) for j in range(dof)]
        dq = [float(numpy.random.random() * 2.0 * pi - pi) for j in range(dof)]
        ddq = [float(numpy.random.random() * 2.0 * pi - pi)
               for j in range(dof)]
        Z[i * dof: i * dof + dof, :] = numpy.matrix(
            regressor_func(q, dq, ddq)).reshape(dof, parm_num)

    R1_diag = numpy.linalg.qr(Z, mode='economic').diagonal().round(round)
    dbi = []
    ddi = []
    for i, e in enumerate(R1_diag):
        if e != 0:
            dbi.append(i)
        else:
            ddi.append(i)
    dbn = len(dbi)

    P = numpy.mat(numpy.eye(parm_num))[:, dbi + ddi]
    Pb = P[:, :dbn]
    Pd = P[:, dbn:]

    Rbd1 = numpy.mat(numpy.linalg.qr(Z * P, mode='r'))
    Rb1 = Rbd1[:dbn, :dbn]
    Rd1 = Rbd1[:dbn, dbn:]

    Kd = numpy.mat((numpy.linalg.inv(Rb1) * Rd1).round(round))

    return Pb, Pd, Kd
