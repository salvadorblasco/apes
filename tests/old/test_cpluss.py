import sys
sys.path.append('../src/')

import pytest as pt
import hypothesis as hp
import hypothesis.strategies as st
import hypothesis.extra.numpy as hen
import numpy as np

from libeq.cpluss import cpluss

@hp.composite
@hp.given(n_points=st.integers(min_value=1, max_value=2000),
          n_equil=st.integers(min_value=1, max_value=500),
          n_species=st.integers(min_value=1, max_value=500))
def generate_cbs(n_points, n_equil, n_species):
    return np.random.rand(n_points, n_species), \
           np.random.rand(n_equil),             \
           np.random.rand(n_equil, n_species)


# @hp.settings(max_examples=500)
# def test_general(n_points, n_equil, n_species):
#     concentration, beta, stoichiometry = generate_cbs(n_points, n_equil, n_species)

@hp.given(n_points=st.integers(min_value=1, max_value=2000),
          n_equil=st.integers(min_value=1, max_value=500),
          n_species=st.integers(min_value=1, max_value=500))
@hp.given(concentration=hen.arrays(np.float, (n_points,n_species), elements=floats(0, 1)),
       beta=hen.arrays(np.float, n_equil, elements=floats(0, 1)),
       stoichiometry=hen.arrays(np.int, (n_equil,n_species), elements=integers())
)
def test_general(concentration, beta, stoichiometry):
    res = cpluss(concentration, beta, stoichiometry)
    assert isinstance(res, np.ndarray)
    assert res.ndim == 2
    assert res.shape == (n_points, n_equil)
    # np.testing.assert_array_almost_equal(res,
