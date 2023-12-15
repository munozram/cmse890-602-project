# Tutorial
*Tutorials are lessons that take the reader by the hand through a series of steps to complete a project of some kind. Tutorials are learning-oriented.*

The code in this project is used to analyze the general structure of a directed graph. Here we show a simple example that uses this tool.

### A weighted directed graph

First, we need a graph to analyze. For this example we will use the following graph given by its adjacency matrix. Note that the weights of edges are included in this matrix.

```
import numpy as np
adjacency_matrix = np.array([[0. , 1. , 0. , 1.6, 0. ],
                            [0. , 0. , 0. , 0. , 0. ],
                            [0. , 2. , 0. , 0. , 0. ],
                            [0. , 0. , 0.8, 0. , 0. ],
                            [0.5, 0. , 0. , 3. , 0. ]])
graph = nx.DiGraph(adjacency_matrix)
```
![Example graph](example_graph.png)

### Compute persistent homology

Use the function as shown below:
```
>>> newfiltration_persistence( graph )
>>> ...
```