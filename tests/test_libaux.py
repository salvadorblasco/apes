#!/usr/bin/python3

import random
import sys
import unittest

import numpy as np

import hypothesis as hp
import hypothesis.strategies as st

sys.path.append('../src/')
import libaux


class LibAuxTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test_assert_array(self):
        a = np.array([1,2,3])
        self.assertRaises(TypeError, libaux.assert_array, 1.0)
        self.assertRaises(TypeError, libaux.assert_array, "a")
        self.assertIsNone(libaux.assert_array(a))

    @hp.given(dim1=st.integers(min_value=1, max_value=10),
              dim2=st.integers(min_value=1, max_value=10))
    def test_assert_array_dim(self, dim1, dim2):
        array_shape = dim1*(2,)
        array = np.empty(array_shape)
        func = libaux.assert_array_dim
        assert array.ndim == dim1
        if dim1 == dim2:
            self.assertIsNone(func(dim2, array))
        else:
            self.assertRaises(ValueError, func, dim2, array)

    @unittest.skip('not implemented')
    def test_assert_shape(self):
        pass

    @unittest.skip('not implemented')
    def test_assert_type(self):
        pass

    @hp.given(list1=st.lists(st.integers()),
              list2=st.lists(st.integers()))
    def test_assert_same_len1(self, list1, list2):
        func = libaux.assert_same_len
        if len(list1) == len(list2):
            self.assertIsNone(func(list1, list2))
        else:
            self.assertRaises(ValueError, func, list1, list2)

    def test_assert_same_shape(self):
        with self.subTest(msg='equal'):
            a = np.random.rand(3,4)
            b = np.random.rand(3,4)
            self.assertTrue(libaux.assert_same_shape(a,b))

        with self.subTest(msg='different'):
            a = np.random.rand(3,4)
            b = np.random.rand(3,2)
            with self.assertRaises(ValueError):
                libaux.assert_same_shape(a,b)

    @unittest.skip('not implemented')
    def test_assert_same_shapeN(self):
        pass

    @unittest.skip('not implemented')
    def test_assert_BPT_consistency(self):
        pass

    # def test_unravel(self):
    #     x = ('a', 'b', 'c', 'd')
    #     test_flags = {
    #         (1,1,1,1): 'abcd',
    #         (0,1,1,1): 'bcd',
    #         (1,1,2,2): 'abc',
    #         (0,0,0,0): ''
    #     }

    #     for flags, result in test_flags.items():
    #         msg = "".join(map(str, flags))
    #         with self.subTest(msg):
    #             retv = "".join(libaux.unravel(x, flags))
    #             self.assertEqual(retv, result)

    # def test_ravel(self):
    #     x = ('a', 'b', 'c', 'd')
    #     y = ('A', 'B', 'C', 'D')
    #     test_flags = {
    #         (1,1,1,1): 'ABCD',
    #         (0,1,1,1): 'aABC',
    #         (0,0,0,0): 'abcd'
    #     }
    #     for flags, result in test_flags.items():
    #         msg = "".join(map(str, flags))
    #         with self.subTest(msg):
    #             retv = "".join(libaux.ravel(x, y, flags))
    #             self.assertEqual(retv, result)

    #     x = (1, 2, 3, 4)
    #     y = (4, 8)
    #     flags = (0, 2, 0, 2)
    #     result = (1, 4, 3, 8)
    #     out = tuple(libaux.ravel(x, y, flags))
    #     self.assertTupleEqual(out, result)

    @unittest.skip('not implemented')
    def test_check_refine_flags(self):
        pass

    @unittest.skip('deprecate')
    def test_count_flags(self):
        test_list = {(0,0,1,1,0,1): (3, 3, {}),
                     (0,0,0,0,0,0): (6, 0, {}),
                     (0,0,1,1,2,2): (2, 2, {2:2})}
        for flags, result in test_list.items():
            msg = "".join(map(str, flags))
            with self.subTest(msg):
                retv = libaux.count_flags(flags)
                self.assertEqual(retv, result)

    @unittest.skip('not implemented')
    def test_build_multiT_titr3(self):
        pass

    @unittest.skip('not implemented')
    def test_build_multiT_titr(self):
        pass

    @unittest.skip('not implemented')
    def test_build_multiT_titr2(self):
        pass

    @unittest.skip('not implemented')
    def test_build_T_titr(self):
        pass

    @unittest.skip('not implemented')
    def test_build_T_titr2(self):
        pass

    @unittest.skip('not implemented')
    def test_build_T(self):
        pass

    def test_extract_groups(self):
        array = np.array([[1, 0, 1],
                          [1, 0, 2],
                          [1, 1, 1],
                          [1, 1, 2],
                          [1, 3, 2],
                          [1, 3, 2]])
        ret = libaux.extract_groups(array, xref=2)
        for a, b in ((0,1), (2,3), (4,5)):
            self.assertEqual(ret[a], ret[b])

        array = np.array([[1, 0, 0, 1],
                          [1, 0, 0, 2],
                          [1, 1, 1, 1],
                          [1, 1, 1, 2],
                          [1, 3, 1, 1],
                          [1, 3, 1, 2]])
        ret = libaux.extract_groups(array, xref=3, yref=2)
        for a, b in ((0,1), (2,3)):
            self.assertEqual(ret[a], ret[b])

    @unittest.skip('not implemented')
    def test_linspace(self):
        pass

    @unittest.skip('not implemented')
    def test_get_lims(self):
        pass

    @unittest.skip('not implemented')
    def test_make_ranges(self):
        pass

    @unittest.skip('not implemented')
    def test_nwe2str(self):
        pass

    @unittest.skip('not implemented')
    def test_nwe2arr(self):
        pass

    @unittest.skip('not implemented')
    def test_numwerr(self):
        pass

    @unittest.skip('not implemented')
    def test_inumwerr(self):
        pass

    def test_percent_distribution(self):
        with np.load('distr_cu2pz32323.npz') as f:
            freec: np.ndarray[float] = f['free_concentration']
            stoich: np.ndarray[int]  = f['stoichiometry']
            analc: np.ndarray[float] = f['analyticalc']

        ref = 1
        x = libaux.percent_distribution(freec, stoich, analc, reference=ref)
        rP = stoich[:, ref].tolist()
        q_plot = [i for i in rP if i != 0]
        q_plot.insert(0, 1)
        st = np.array(q_plot)
        y = np.sum(x, axis=1)
        # y = np.sum(x*st[None,:],axis=1)
        np.testing.assert_array_almost_equal(np.full_like(y, 100.0), y)


if __name__ == '__main__':
    unittest.main()
