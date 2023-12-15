import pytest
import numpy as np
import networkx as nx
from PHN_Directed import all_fvalues

# Cyclic graph
A1 = np.diag(np.ones((5,)),1)
A1[-1,0] = 1
D1 = nx.DiGraph(A1)
result1 = np.array([0., 0., 0., 0., 0., 0., 1., 2., 3., 2., 1., 1., 2., 3., 2., 1., 2.,
       3., 1., 2., 1., 2., 3., 3., 2., 3., 4., 3., 3., 3., 2., 2., 3., 3.,
       3., 4., 3., 2., 3., 3., 2., 3., 4., 3., 4., 4., 3., 4., 4., 4., 3.,
       3., 4., 4., 4., 3., 4., 4., 4., 4., 4., 4., 5.])

# Almost cyclic graph
A2 = np.diag(np.ones((5,)),1)
A2[-1,0] = 6
A2[0,-1] = 1
D2 = nx.DiGraph(A2)
result2 = np.array([0., 0., 0., 0., 0., 0., 1., 2., 3., 4., 1., 1., 2., 3., 4., 1., 2.,
       3., 1., 2., 1., 2., 3., 4., 5., 3., 4., 5., 4., 5., 5., 2., 3., 4.,
       3., 4., 4., 2., 3., 3., 2., 3., 4., 5., 4., 5., 5., 4., 5., 5., 5.,
       3., 4., 4., 4., 3., 4., 5., 5., 5., 5., 4., 5.])

# Random graph
A3 = np.array([[0. , 1. , 0. , 1.6, 0. ],
               [0. , 0. , 0. , 0. , 0. ],
               [0. , 2. , 0. , 0. , 0. ],
               [0. , 0. , 0.8, 0. , 0. ],
               [0.5, 0. , 0. , 3. , 0. ]])
D3 = nx.DiGraph(A3)
result3 = np.array([0. , 0. , 0. , 0. , 0. , 1. , 2.4, 1.6, 0.5, 2. , 2.8, 1.5, 0.8,
       3.8, 3. , 4.4, 4.4, 1.5, 2.4, 2.9, 2.1, 2.8, 5.8, 5.8, 3.8, 4.4,
       4.9, 4.9, 2.9, 5.8, 4.9])

@pytest.mark.parametrize("D,result", [(D1,result1),(D2,result2),(D3,result3)])
def test_filtration(D,result):
    _, fvals = all_fvalues(D, max_size=7)
    np.testing.assert_array_almost_equal(fvals, result)
