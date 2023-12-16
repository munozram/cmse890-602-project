# Test for the function `newfiltration_persistence`
# Should only be run locally where the python package Dionysus is installed

import pytest
from PHN_Directed import *

# Cyclic graph
A1 = np.diag(np.ones((5,)), 1)
A1[-1, 0] = 1
D1 = nx.DiGraph(A1)
result1 = np.array([[0., 0., 0., 0., 0., 0., 1., 3., 3., 4.],
                    [np.inf, 1., 1., 1., 1., 1., 3., 4., 4., 5.],
                    [0., 0., 0., 0., 0., 0., 1., 2., 2., 4.]])

# Almost cyclic graph
A2 = np.diag(np.ones((5,)), 1)
A2[-1, 0] = 6
A2[0, -1] = 1
D2 = nx.DiGraph(A2)
result2 = np.array([[0., 0., 0., 0., 0., 0., 1.],
                    [np.inf, 1., 1., 1., 1., 1., 5.],
                    [0., 0., 0., 0., 0., 0., 1.]])

# Random graph
A3 = np.array([[0., 1., 0., 1.6, 0.],
               [0., 0., 0., 0., 0.],
               [0., 2., 0., 0., 0.],
               [0., 0., 0.8, 0., 0.],
               [0.5, 0., 0., 3., 0.]])
D3 = nx.DiGraph(A3)
result3 = np.array([[0., 0., 0., 0., 0., 2., 2.1, 2.9, 2.9, 4.9, 4.9, 4.9],
                    [np.inf, 1., 1.6, 0.8, 0.5, 4.4, 3., 3.8, 3.8, 5.8, 5.8, 5.8],
                    [0., 0., 0., 0., 0., 1., 2., 2., 3., 3., 3., 4.]])


@pytest.mark.parametrize("D,result",
                         [(D1, result1), (D2, result2), (D3, result3)])
def test_filtration(D, result):
    dgms = newfiltration_persistence(D, max_dim=4)
    births, deaths, dims = [], [], []
    for i, dgm in enumerate(dgms):
        births = births + [pt.birth for pt in dgm]
        deaths = deaths + [pt.death for pt in dgm]
        dims = dims + [i for _ in dgm]
    this_result = np.array([births, deaths, dims])
    np.testing.assert_array_almost_equal(this_result, result)
