#!/usr/bin/python3

import unittest

import sys
import numpy as np

# import hexaprotic

sys.path.append('../src')
# pylint: disable=import-error


class Test0General(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test0_engine_exists(self):
        import libeq
        self.assertTrue('LIBEQ_ENGINE' in dir(libeq))

    def test1_engine_valid(self):
        import libeq
        self.assertIn(libeq.LIBEQ_ENGINE, ('python', 'fortran'))


class Test1Cpluss(unittest.TestCase):
    # TODO test logc option
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test_python_cpluss(self):
        from libeq.cpluss import cpluss
        beta, stoichiometry, analyticalc = self.__hardest()
        n_equil, n_species = stoichiometry.shape
        n_points = analyticalc.shape[0]
        conc = np.random.rand(n_points, n_species)
        econc = cpluss(conc, beta, stoichiometry)
        self.__common_cpluss(beta, stoichiometry, analyticalc, conc, econc)

    def test_python_cpluss_masked(self):
        from libeq.cpluss import cpluss
        beta, stoichiometry, analyticalc = self.__hardest()
        n_equil, n_species = stoichiometry.shape
        n_points = analyticalc.shape[0]
        import random
        _mask = np.column_stack(n_species*[random.choice([True, False])
                                           for _ in range(n_points)])
        conc = np.ma.array(np.random.rand(n_points, n_species), mask=_mask)
        econc = cpluss(conc, beta, stoichiometry)
        self.__common_cpluss(beta, stoichiometry, analyticalc, conc, econc)

    def test_fortran_cpluss(self):
        import libeq.mfeq
        beta, stoichiometry, analyticalc = self.__hardest()
        n_equil, n_species = stoichiometry.shape
        n_points = analyticalc.shape[0]

        f_beta = np.asfortranarray(beta)
        f_stoichiometry = np.asfortranarray(stoichiometry)
        f_cextd = np.empty((n_points, n_equil+n_species), order='F')
        f_cextd[:, :n_species] = np.random.rand(n_points, n_species)
        libeq.mfeq.cplussn(f_cextd, f_beta, f_stoichiometry, n_points, n_equil,
                           n_species)
        self.__common_cpluss(beta, stoichiometry, analyticalc,
                             f_cextd[:, :n_species], f_cextd[:, n_species:])

    def __common_cpluss(self, beta, stoichiometry, analyticalc, conc, econc):
        def prod(numbers):
            from functools import reduce
            return reduce(lambda a, b: a*b, numbers)

        n_equil, n_species = stoichiometry.shape
        n_points = analyticalc.shape[0]
        f_tests = 0.5   # test 50% of points
        n_tests = int(f_tests*n_points*n_equil)
        for n in range(n_tests):
            r = np.random.randint(0, n_points)
            c = np.random.randint(0, n_equil)
            _ = prod(conc[r, j]**stoichiometry[c, j] for j in range(n_species))
            if not np.ma.is_masked(_):
                comp = econc[r, c] / (beta[c]*_)
                self.assertAlmostEqual(comp, 1.0)

    def __hardest(self):
        import libaux
        labels = ('L', 'Cu', 'Zn', 'Im', 'H')
        S = len(labels)
        log10_beta = np.array((10.97, 20.52, 29.12, 36.61, 43.73, 48.72,
                               7.09, 4.2, 7.73, 10.6, 13.0, 44.77, 40.08,
                               35.43, 27.87, 18.34, 30.03, 22.26, 39.2, 33.56,
                               25.2, 17.08, 33.66, 26.69, 20.04, 10.76, 15.69,
                               8.1, -1.8, 26.03, 18.36, 11.22, 2.47, 23.22,
                               15.87, 5.28, 5.28, -13.73))
        E = len(log10_beta)
        stoichiometry = np.array((1, 0, 0, 0, 1,
                                  1, 0, 0, 0, 2,
                                  1, 0, 0, 0, 3,
                                  1, 0, 0, 0, 4,
                                  1, 0, 0, 0, 5,
                                  1, 0, 0, 0, 6,
                                  0, 0, 0, 1, 1,
                                  0, 1, 0, 0, 1,
                                  0, 1, 0, 0, 2,
                                  0, 1, 0, 0, 3,
                                  0, 1, 0, 0, 4,
                                  1, 1, 0, 0, 4,
                                  1, 1, 0, 0, 3,
                                  1, 1, 0, 0, 2,
                                  1, 1, 0, 0, 1,
                                  1, 1, 0, 0, 0,
                                  1, 2, 0, 0, 0,
                                  1, 2, 0, 0, -1,
                                  1, 2, 0, 1, 1,
                                  1, 2, 0, 1, 0,
                                  1, 2, 0, 1, -1,
                                  1, 2, 0, 1, -2,
                                  1, 0, 1, 0, 3,
                                  1, 0, 1, 0, 2,
                                  1, 0, 1, 0, 1,
                                  1, 0, 1, 0, 0,
                                  1, 0, 2, 0, 0,
                                  1, 0, 2, 0, -1,
                                  1, 0, 2, 0, -2,
                                  1, 0, 2, 1, 1,
                                  1, 0, 2, 1, 0,
                                  1, 0, 2, 1, -1,
                                  1, 0, 2, 1, -2,
                                  1, 1, 1, 0, 0,
                                  1, 1, 1, 0, -1,
                                  1, 1, 1, 0, -2,
                                  1, 1, 1, 1, -1,
                                  0, 0, 0, 0, -1)).reshape((E, S))
        t0 = 0.01, 0.01, 0.01, 0.01, 0.10
        buret = 0.0, 0.0, 0.0, 0.0, -0.1
        v0 = 10.000000
        v1 = 11.000000
        N = 10
        incrV = (v1-v0)/N
        beta = 10**log10_beta
        T = np.array(tuple(libaux.build_T_titr(t0, buret, v0, incrV, N)))
        return beta, stoichiometry, T


class Test2Jacobian(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.logB = np.array([10.13, 19.53, 27.80, 34.82, -13.73])
        self.B = 10**self.logB
        self.P = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [0, -1]])
        with np.load('pytrenal.npz') as d:
            self.real_C: np.ndarray[float] = d['C']
            self.real_T: np.ndarray[float] = d['T']
        self.E, self.S = self.P.shape

    def test_python_jacobian(self):
        from libeq.jacobian import jacobian
        J = jacobian(self.real_C, self.P)
        n = self.real_C.shape[0]
        self.assertEqual(J.shape, (n, self.S, self.S))
        for m in range(50):
            i, jvalue = randompick(J)
            jcalc = sum((self.P[k, i[1]]*self.P[k, i[2]]
                         *self.real_C[i[0], k+self.S]
                         for k in range(self.E)))
            jcalc /= self.real_C[i[0], i[2]]
            jcalc += 1 if i[1] == i[2] else 0
            self.assertAlmostEqual(jcalc, jvalue)

    def test_fortran_jacobian(self):
        import libeq.mfeq
        c = np.empty(self.E + self.S, order='F')
        jin = np.empty((self.S, self.S), order='F')
        p = np.asfortranarray(self.P)
        for _c in self.real_C:
            c[:] = _c
            libeq.mfeq.jacobian(jin, c, p, self.E, self.S)
            for i in range(self.S):
                for j in range(self.S):
                    jout = sum((self.P[k, i]*self.P[k, j]*c[k+self.S]
                                for k in range(self.E)))
                    jout /= c[j]
                    jout += 1 if i == j else 0
                    self.assertAlmostEqual(jin[i, j], jout)

    def test_amatrix(self):
        import libeq.jacobian
        import hexaprotic
        amatrix = libeq.jacobian.amatrix(hexaprotic.free_concentration, hexaprotic.stoichx)
        np.testing.assert_array_almost_equal(amatrix, hexaprotic.matrix_a, decimal=4)

    def test_amatrix_lmh(self):
        import libeq.jacobian
        import data_lmh

        with self.subTest(dataset='t1'):
            amatrix_t1 = libeq.jacobian.amatrix(data_lmh.t1_freeconc, data_lmh.stoichx)
            np.testing.assert_array_almost_equal(amatrix_t1, data_lmh.t1_amatrix)

        with self.subTest(dataset='t2'):
            amatrix_t2 = libeq.jacobian.amatrix(data_lmh.t2_freeconc, data_lmh.stoichx)
            np.testing.assert_array_almost_equal(amatrix_t2, data_lmh.t2_amatrix)

    def test_dlogcdlogbeta(self):
        import libeq.jacobian
        import hexaprotic
        tested = libeq.jacobian.dlogcdlogbeta(hexaprotic.matrix_a,
                                              hexaprotic.free_concentration,
                                              hexaprotic.stoich)
        np.testing.assert_array_almost_equal(tested, hexaprotic.dlogc_dlogbeta, decimal=2)

    def test_dlogcdt(self):
        import libeq.jacobian
        import hexaprotic
        tested = libeq.jacobian.dlogcdt(hexaprotic.matrix_a,
                                        hexaprotic.titre,
                                        hexaprotic.v0)
        np.testing.assert_allclose(tested, hexaprotic.dlogc_dt, rtol=1e-2)

    def test_dlogcdb(self):
        import libeq.jacobian
        import hexaprotic
        tested = libeq.jacobian.dlogcdb(hexaprotic.matrix_a,
                                        hexaprotic.titre,
                                        hexaprotic.v0)
        np.testing.assert_allclose(tested, hexaprotic.dlogc_db, rtol=1e-2)



class Test3Fobj(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.logB = np.array([10.13, 19.53, 27.80, 34.82, -13.73])
        self.B = 10**self.logB
        self.P = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [0, -1]])
        with np.load('pytrenal.npz') as d:
            self.real_C = d['C']
            self.real_T = d['T']
        self.E, self.S = self.P.shape

    def test_python_fobj(self):
        from libeq.fobj import fobj
        f = fobj(self.real_C, self.P, self.real_T)
        np.testing.assert_array_less(np.abs(f), np.full_like(f, 1e-5))

    def test_fortran_fobj(self):
        import libeq.mfeq
        f_conc = np.asfortranarray(self.real_C)
        f_stoi = np.asfortranarray(self.P)
        f_anac = np.asfortranarray(self.real_T)
        F = np.empty(self.S, order='F')
        for i in range(self.real_C.shape[0]):
            libeq.mfeq.fobj(F, f_conc[i], f_stoi, f_anac[i], self.E, self.S)
            np.testing.assert_array_less(np.abs(F), np.full_like(F, 1e-5))


class Test4Consol0(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.logB = np.array([10.13, 19.53, 27.80, 34.82, -13.73])
        self.B = 10**self.logB
        self.P = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [0, -1]])
        with np.load('pytrenal.npz') as d:
            self.real_C = d['C']
            self.real_T = d['T']
        self.E, self.S = self.P.shape

    @unittest.skip('broken')
    def test_python_consol(self):
        import libeq
        libeq.LIBEQ_ENGINE = 'python'
        import libeq.consol
        for T in map(np.array, self.real_T):
            # breakpoint()
            x0 = libeq.consol.initial_guess(self.B, self.P, self.real_T)
            calc_C = libeq.consol.consol(self.B, self.P, self.real_T, x0) #[:,:self.S])
            np.testing.assert_allclose(self.real_C, calc_C)

    @unittest.skip('broken')
    def test_fortran_consol(self):
        import libeq
        libeq.LIBEQ_ENGINE = 'fortran'
        import libeq.consol
        for T in map(np.array, self.real_T):
            calc_C = libeq.consol.consol(self.B, self.P, self.real_T,
                                         self.real_C[:, :self.S])
            np.testing.assert_allclose(self.real_C, calc_C)


class Test4Initial(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test_simu(self):
        import libeq.consol
        filename = 'distr_phos.npz'
        with np.load(filename) as f:
            realc = f['free_concentration']
            beta = 10**f['beta']
            analc = f['analyticalc']
            stoich = f['stoichiometry']

        # n_points = len(realc)
        for engine in ('python', 'fortran'):
            libeq.LIBEQ_ENGINE = engine
            with self.subTest('engine = ' + engine):
                new_x, beta_prime, stoich_new, analc_new = \
                    libeq.consol.freeze_concentration(beta, stoich, analc,
                                                      reference=1)
                conc = libeq.consol.initial_guess(beta_prime, stoich_new,
                                                  analc_new)
                np.testing.assert_array_almost_equal(np.delete(realc, 1, axis=1), conc, decimal=2)
        

class Test6Panic(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def feed_files(self):
        import glob
        yield from glob.iglob('consol_panic_*.npz')

    def run_initial_guess(self, filename):
        import numpy as np
        from libeq.consol import initial_guess
        with np.load(filename) as f:
            beta = f['beta']
            analytc = f['analytc']
            stoich = f['stoichiometry']
        initial_guess(beta, stoich, analytc, panic=False)

    @unittest.skip('broken')
    def test_panic_files_fortran(self):
        # import glob
        # import numpy as np
        # from libeq.consol import initial_guess
        import libeq
        libeq.LIBEQ_ENGINE = 'fortran'
        for filename in self.feed_files():
            with self.subTest():
                # with np.load(filename) as f:
                #     beta = f['beta']
                #     analytc = f['analytc']
                #     stoich = f['stoichiometry']
                # initial_guess(beta, stoich, analytc, panic=False)
                self.run_initial_guess(filename)

    @unittest.skip('broken')
    def test_panic_files_python(self):
        # import glob
        # import numpy as np
        # from libeq.consol import initial_guess
        import libeq
        libeq.LIBEQ_ENGINE = 'python'
        # for filename in glob.iglob('consol_panic_*.npz'):
        for filename in self.feed_files():
            with self.subTest():
                # with np.load(filename) as f:
                #     beta = f['beta']
                #     analytc = f['analytc']
                #     stoich = f['stoichiometry']
                # initial_guess(beta, stoich, analytc, panic=False)
                self.run_initial_guess(filename)


class Test7Solid(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        pass

    @unittest.skip('incomplete')
    def test_g(self):
        # TODO. Complete.
        log_beta1 = np.array([-9.0, -18.7, -27.0, -33.0, -13.73])
        stoich1 = np.array([[1,-1],[1,-2],[1,-3],[1,-4],[0,-1]])

        log_beta2 = np.array([-33.5])
        stoich2 = np.array([[1,-3]])


# def random_mask(array, axis=0):
#     #TODO
#     import random
#     _mask = np.column_stack(n_species*[random.choice([True, False])
#                                            for _ in range(n_points)])
#     return np.ma.array(array, mask=_mask)


def randompick(array):
    pick = tuple(np.random.randint(0, dim) for dim in array.shape)
    return pick, array[pick]


if __name__ == '__main__':
    unittest.main()
