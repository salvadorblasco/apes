#!/usr/bin/python3

import unittest

import numpy as np

from context import libeq, libaux, libemf


def ccalc(T, B, Kw=1e-14):
    # import pudb
    # pudb.set_trace()
    C = np.zeros((T.shape[0], len(B)+3))
    H = np.fromiter(hcalc(T, B, Kw), dtype=np.float)
    L = T[:, 0]/(1+np.sum(B[np.newaxis, :]*H[:, np.newaxis]**np.arange(1, 1+len(B)), axis=1))
    C[:, 0] = L
    C[:, 1] = H
    for i, b in enumerate(B, start=2):
        C[:, i] = b*L*H**(i-1)
    C[:, -1] = Kw/H
    return C


def hcalc(T, B, Kw):
    # degree = len(B) + 2
    coeffs = (len(B) + 3)*[0.0]
    for t in T:
        for i in range(len(B)):
            coeffs[i+2] += B[i]
            coeffs[i+1] += t[0]*i*B[i]
            coeffs[i+1] -= t[1]*B[i]
            coeffs[i] -= Kw*B[i]
        coeffs[2] += 1
        coeffs[1] -= t[0]
        coeffs[0] -= Kw
        r = np.roots(tuple(reversed(coeffs)))
        h = filter_h(r)
        yield h


def filter_h(roots):
    be_real = [x for x in roots if isinstance(x, (float, int)) or x.imag == 0]
    be_positive = [x for x in be_real if x >= 0.0]
    if len(be_positive) == 1:
        return be_positive[0]
    else:
        return min(be_positive)


class GenericProtonation(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        T0 = np.array([0.001, 0.002])
        V0 = 20.0
        buret = np.array([0.0, -0.1])
        N = 20
        V = np.linspace(0.0, 1.2, N)
        self.T = np.array(libaux.build_T_titr2(T0, buret, V0, V))

    def test_monoprotic(self):
        P = np.array([[1, 1], [0, -1]])
        B = 10**np.array([9.75, -13.73])
        real_C = ccalc(self.T, B[:-1], Kw=B[-1])
        calc_C = libeq.consol(B, P, self.T)
        #import pudb
        #pudb.set_trace()
        # from matplotlib import pyplot
        # pyplot.plot(real_C, color='red')
        # pyplot.plot(calc_C, color='blue')
        # pyplot.figure()
        # #pyplot.plot(np.log10(real_C[:,2]/(real_C[:,0]*real_C[:,1])), color='red')
        # pyplot.plot(np.log10(real_C[:,3]*real_C[:,1]), color='red')
        # #pyplot.plot(np.log10(calc_C[:,2]/(calc_C[:,0]*calc_C[:,1])), color='blue')
        # pyplot.plot(np.log10(calc_C[:,3]*calc_C[:,1]), color='blue')
        # pyplot.show()
        print(libeq.fobj(real_C, P, self.T))
        print(libeq.fobj(calc_C, P, self.T))
        print(np.sum(libeq.fobj(calc_C, P, self.T)))
        np.testing.assert_allclose(real_C, calc_C)

    def test_zeroprotic(self):
        pKw = -13.73
        P = np.array([[-1]])
        B = 10**np.array([pKw])
        T = self.T[:, 1, np.newaxis]
        calc_C = libeq.consol(B, P, T)
        aux = np.prod(calc_C, axis=1)
        np.testing.assert_allclose(np.log10(aux),
                                   pKw*np.ones_like(aux))


class TestAcetic(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def acetic_realC(self, T, K=None):
        if K is None:
            K = self.B[0]
        Kw = self.B[1]
        N = T.shape[0]
        C = np.empty((N, 4))
        a = K
        b = 1+K*(T[:, 0] - T[:, 1])
        c = -T[:, 1]-K*Kw
        d = -Kw
        for n in range(N):
            r = np.roots([a, b[n], c[n], d])
            h = r[r[np.isreal(r)] >= 0][0]
            assert h >= 0
            C[n, 1] = h
            C[n, 0] = T[n, 0] / (1+K*h)
            C[n, 2] = T[n, 0] - C[n, 0]
            C[n, 3] = Kw/h
        return C

    def setUp(self):
        self.logB = np.array([4.76, -13.73])
        self.B = 10**self.logB
        self.P = np.array([[1, 1], [0, -1]])
        self.E, self.S = self.P.shape

    def test_acetic(self):
        # simple equilibrium
        # Kh^3 + (1+K(T_L - T_H ))h^2-(T_H+KK_w )h - K_w = 0
        T0 = np.array([0.1, 0.1])
        V0 = 20.0
        buret = np.array([0.0, -0.1])
        N = 20
        V = np.linspace(0.0, 1.0, N)
        T = np.array(libaux.build_T_titr2(T0, buret, V0, V))
        real_C = self.acetic_realC(T)
        calc_C = libeq.consol(self.B, self.P, T)
        np.testing.assert_allclose(real_C, calc_C)

    def test_dcdB_acetic(self):
        # simple equilibrium
        # Kh^3 + (1+K(T_L - T_H ))h^2-(T_H+KK_w )h - K_w = 0
        T0 = np.array([0.1, 0.1])
        V0 = 20.0
        buret = np.array([0.0, -0.1])
        N = 20
        V = np.linspace(0.0, 1.0, N)
        T = np.array(libaux.build_T_titr2(T0, buret, V0, V))
        C = self.acetic_realC(T)

        calc_dcdB = libeq.extdd_dcdb(C, self.P)
        real_dcdB = np.empty_like(calc_dcdB)
        for n, conc in enumerate(C):
            l, h, c, z = conc
            dhdK = -c*l / (c*h + c*z + l*h + l*c + l*z)
            dhdKw = z*(c+l) / (c*h + c*z + l*h + l*c + l*z)
            dldK = -c/(l+c)*(1+dhdK)
            dldKw = -c/(l+c)*dhdKw
            dcdK = 1 + dhdK + dldK
            dcdKw = dhdKw+dldKw
            dzdK = -dhdK
            dzdKw = 1-dhdKw

            real_dcdB[n] = np.array([[dldK,  dhdK,  dcdK, dzdK],
                                     [dldKw, dhdKw, dcdKw, dzdKw]]).T

        np.testing.assert_almost_equal(real_dcdB, calc_dcdB)

    def test_jacK_acetic(self):
        """simple equilibrium
        Kh^3 + (1+K(T_L - T_H ))h^2-(T_H+KK_w )h - K_w = 0
        dh/dK = {-h^3 - (T_L-T_H)h^2 + K_w h} /
                {3Kh^2 + 2[1+K(T_L - T_H)]h - T_H+KK_w}
        J = {d log h}/{d log K} = K/h dh/dK =
          = { -h^2 -( T_L-T_H )h + K_w  } /
            { 3h^2 + 2[1/K+T_L-T_H]h - T_H/K - K_w }
        """
        T0 = np.array([0.1, 0.2])
        V0 = 20.0
        buret = np.array([0.0, -0.1])
        N = 20
        V = np.linspace(0.0, 2.5, N)
        T = np.array(libaux.build_T_titr2(T0, buret, V0, V))
        real_C = np.empty((N, 4))
        # E0 = 380.0  # mV
        # d = -Kw

        for lK in (4.76, 4.50, 4.76, 4.95):
            real_C = self.acetic_realC(T, K=10**lK)

            H = real_C[:, 1]
            real_J = -real_C[:, 2]*real_C[:, 0] / (
                real_C[:, 2]*H +
                real_C[:, 0]*H +
                real_C[:, 2]*real_C[:, 0] +
                real_C[:, 2]*real_C[:, 3] +
                real_C[:, 0]*real_C[:, 3])
            assert len(real_J) == N

            calc_C = libeq.consol(self.B, self.P, T, x0=real_C[:, :2])
            np.testing.assert_allclose(real_C, calc_C)

            calc_J = libemf.emf_jac1(self.P, calc_C, 1, (0,)).reshape(N)
            np.testing.assert_allclose(real_J, calc_J)


class TestHPytren(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test_jacK_hpytren(self):
        """# simple equilibrium
        # B3h^5 + ah^4 + bh^3 + ch^2 + dh -Kw = 0
        #   a = B3(3TL-TH) + B2
        #   b = B2(2TL-TH) - B3.Kw + B1
        #   c = 1+B1(TL-TH) - B2.Kw
        #   d = TH-B1.Kw
        # dh/dB1 = (-h^3-(TL-TH)h^2+Kw.h) / D
        # dh/dB2 = (-h^4-(2.TL-TH)h^3+Kw.h^2) / D
        # dh/dB3 = (-h^5-(3.TL-TH)h^4+Kw.h^3) / D
        #   D = 5.B3.h^4 + 4a.h^3 + 3b.h^2 + 2c.h + d"""
        lKw = -13.73
        Kw = 10.0**(lKw)
        P = np.array([[1, 1], [1, 2], [1, 3], [0, -1]])
        E, S = P.shape
        T0 = np.array([0.1, 0.5])
        V0 = 20.0
        buret = np.array([0.0, -0.1])
        N = 20
        V = np.linspace(0.0, 5.0, N)
        T = np.array(libaux.build_T_titr2(T0, buret, V0, V))
        real_C = np.empty((N, E+S))
        # E0 = 382.55  # mV
        logB = np.array([10.0, 20.0, 27.0, lKw])
        B = np.power(10.0, logB)

        a = B[2]*(3*T[:, 0]-T[:, 1]) + B[1]
        b = B[1]*(2*T[:, 0] - T[:, 1]) - B[2]*Kw + B[0]
        c = 1 + B[0]*(T[:, 0]-T[:, 1]) - B[1]*Kw
        d = -T[:, 1]-B[0]*Kw
        for n in range(N):
            r = np.roots([B[2], a[n], b[n], c[n], d[n], -Kw])
            h = r[np.logical_and(np.isreal(r), r.real >= 0)][0].real
            assert isinstance(h, float) and h >= 0, h
            real_C[n, 1] = h
            l = T[n, 0] / (1+B[0]*h+B[1]*h**2+B[2]*h**3)
            real_C[n, 0] = l
            real_C[n, 2] = B[0]*l*h
            real_C[n, 3] = B[1]*l*h**2
            real_C[n, 4] = B[2]*l*h**3
            real_C[n, 5] = Kw/h

        H = real_C[:, 1]
        D = 5*B[2]*H**4 + 4*a*H**3 + 3*b*H**2 + 2*c*H + d
        real_J = np.empty((N, 3))
        real_J[:, 0] = B[0]*(-H**2-(T[:, 0]-T[:, 1])*H+Kw)/D
        real_J[:, 1] = B[1]*(-H**3-(2*T[:, 0]-T[:, 1])*H**2+Kw*H)/D
        real_J[:, 2] = B[2]*(-H**4-(3*T[:, 0]-T[:, 1])*H**3+Kw*H**2)/D

        calc_C = libeq.consol(B, P, T, max_iterations=100, x0=real_C[:, :S])
        # calc_C = real_C

        from libemf import emf_jac1
        calc_J = emf_jac1(P, calc_C, 1, (0, 1, 2))
        assert len(calc_J) == N

        np.testing.assert_almost_equal(real_J, calc_J, decimal=3)


class TestHPytrenAl(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestHPytrenAl, self).__init__(*args, **kwargs)

    def test_cpluss(self):
        C = self.calcC(self.real_T)
        np.testing.assert_allclose(
            C[:, self.S:],
            libeq.cplusS(C[:, :self.S], self.B, self.P))

    def test_consol(self):
        mT = libaux.build_multiT_titr2(self.T0, self.buret, self.V0, self.V)
        for T in map(np.array, mT):
            calc_C = libeq.consol(self.B, self.P, T)
            real_C = self.calcC(T)
            np.testing.assert_allclose(real_C, calc_C)

    def test_jacobian(self):
        T = self.real_T
        N = T.shape[0]
        E, S = self.P.shape
        real_C = self.calcC(T)
        Kw = self.Kw
        B = self.B[:-1]
        a = B[3]*(4*T[:, 0] - T[:, 1]) + B[2]
        b = B[2]*(3*T[:, 0] - T[:, 1]) - B[3]*Kw + B[1]
        c = B[1]*(2*T[:, 0] - T[:, 1]) - B[2]*Kw + B[0]
        d = 1 + B[0]*(T[:, 0] - T[:, 1]) - B[1]*Kw
        e = -T[:, 1]-B[0]*Kw

        H = real_C[:, 1]
        D = 6*B[3]*H**5 + 5*a*H**4 + 4*b*H**3 + 3*c*H**2 + 2*d*H + e
        real_J = np.empty((N, 4))
        real_J[:, 0] = B[0]*(-H**2-(T[:, 0]-T[:, 1])*H+Kw)/D
        real_J[:, 1] = B[1]*(-H**3-(2*T[:, 0]-T[:, 1])*H**2+Kw*H)/D
        real_J[:, 2] = B[2]*(-H**4-(3*T[:, 0]-T[:, 1])*H**3+Kw*H**2)/D
        real_J[:, 3] = B[3]*(-H**5-(3*T[:, 0]-T[:, 1])*H**4+Kw*H**3)/D

        np.testing.assert_almost_equal(real_J,
                                       libeq.jacobian(real_C, self.P, log=True))

    def calcC(self, T, B=None):
        if B is None:
            B = self.B
        Kw = B[-1]
        N = T.shape[0]
        a = B[3]*(4*T[:, 0] - T[:, 1]) + B[2]
        b = B[2]*(3*T[:, 0] - T[:, 1]) - B[3]*Kw + B[1]
        c = B[1]*(2*T[:, 0] - T[:, 1]) - B[2]*Kw + B[0]
        d = 1 + B[0]*(T[:, 0] - T[:, 1]) - B[1]*Kw
        e = -T[:, 1]-B[0]*Kw
        C = np.empty((N, self.E+self.S))
        for n in range(N):
            r = np.roots([B[3], a[n], b[n], c[n], d[n], e[n], -Kw])
            h = r[np.logical_and(np.isreal(r), r.real >= 0)][0].real
            assert isinstance(h, float) and h >= 0, h
            C[n, 1] = h
            l = T[n, 0] / (1+B[0]*h+B[1]*h**2+B[2]*h**3)
            C[n, 0] = l
            C[n, 2] = B[0]*l*h
            C[n, 3] = B[1]*l*h**2
            C[n, 4] = B[2]*l*h**3
            C[n, 5] = B[3]*l*h**4
            C[n, 6] = Kw/h
        return C

    def setUp(self):
        self.T0 = [np.array([0.02521, 0.12606]),
                   np.array([0.01825, 0.09125]),
                   np.array([0.02039, 0.10196])]
        self.V0 = 3*[30.0]
        self.buret = 3*[np.array([0.0, -0.1])]
        self.logB = np.array([10.13, 19.53, 27.80, 34.82, -13.73])
        self.B = 10**self.logB
        self.Kw = self.B[-1]
        self.Bkeys = np.array([0, 0, 0, 0, 1])
        self.P = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [0, -1]])
        self.error_V0 = 3*[0.0030]
        self.E0 = [400.99, 400.66, 397.32]
        self.error_E0 = 3*[0.3000]
        self.V = [np.arange(0.0, 1.301, 0.01),
                  np.arange(0.0, 1.101, 0.01),
                  np.arange(0.0, 1.201, 0.02)]
        with np.load('pytrenal.npz') as d:
            self.real_C = d['C']
            self.real_T = d['T']
            self.emf = [d['emf1'], d['emf2'], d['emf3']]
        self.E, self.S = self.P.shape


if __name__ == '__main__':
    unittest.main()
