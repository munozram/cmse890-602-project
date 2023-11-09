import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import copy

import dionysus as dio

#from PHN_Dowker import *



def pairs_dictionary(D):
    dictionary = {}
    nodes = sorted(list(D.nodes))

    for i,u in enumerate(nodes[:-1]):
        for v in nodes[i+1:]:
            path = None
            length = np.inf
            #print(u,v)
            try:
                path = nx.shortest_path(D,u,v, weight='weight')
                #length = len(path)
                length = nx.path_weight(D, path, weight='weight')
                #print(path)
            except nx.NetworkXNoPath: pass
            try:
                path2 = nx.shortest_path(D,v,u, weight='weight')
                #if len(path2) < length:
                if nx.path_weight(D, path2, weight='weight') < length:
                    path = path2
                #print(path2)
            except nx.NetworkXNoPath: pass
            dictionary[(u,v)] = path
    return dictionary



def find_f_path(D, sigma, dictionary):
    """
    This function finds the shortest path, in a digraph D, that contains all vertices in subset sigma,
    so that the filtration value of sigma is len(path)-1.

    Args:
        D (networkx.DiGraph): directed graph
        sigma (list or tuple): subset of Vert(D)
        dictionary (dictionary): dictionary where shortest paths are stored successively; 
            should already contain at least all shortest paths of all pairs in D.

    Returns:
        [None]: None if there is no path that contains all the vertices
        or
        [list]: shortest path in D containing sigma
    """

    sigma = sorted(sigma)
    try:
        return dictionary[tuple(sigma)]
    except KeyError:
        pass
    length = np.inf

    for i in range(len(sigma)):
        ## remove i from sigma
        subset = copy.copy(sigma)
        subset.remove(sigma[i])
        #print('subset is', subset)

        ## find min path for subset
        path1 = find_f_path(D, subset, dictionary)

        ## if there is no path
        if path1 == None:
            #print('There is no path that contains all the vertices', sigma)
            return None
        else:
            #print('1--', path1)

            ## if i is already in the path, no need to connect i
            if sigma[i] not in path1:
                ## find shortest path to connect i either at the start or at the end
                option = 0
                length2 = np.inf
                try:
                    path2 = nx.shortest_path(D, sigma[i],path1[0], weight='weight')
                    option = 1
                    #length2 = len(path2)
                    length2 = nx.path_weight(D, path2, weight='weight')
                    #print('option 1', path2)
                except nx.NetworkXNoPath: pass
                try:
                    path3 = nx.shortest_path(D, path1[-1],sigma[i], weight='weight')
                    length3 = nx.path_weight(D, path3, weight='weight')
                    #if len(path3) < length2:
                    if length3 < length2:
                        option = 2
                        #path = path2
                    #print('option 2', path3)
                except nx.NetworkXNoPath: pass

                ## concatenate shortest option
                if option == 0:
                    #print('no option')
                    continue
                if option == 1:
                    path1 = path2[:-1] + path1
                if option == 2:
                    path1 = path1 + path3[1:]
            #print('path is', path1)

            ## overwrite path if this is better than the last one
            #if len(path1) < length:
            if nx.path_weight(D, path1, weight='weight') < length:
                #print('hey')
                path = path1
                #length = len(path)
                length = nx.path_weight(D, path, weight='weight')
                if length == len(sigma): break
    #print('completed:', path)
    ## save data in dictionary
    dictionary[tuple(sigma)] = path
    return path


def all_fvalues(D, dictionary, max_size=3):
    """This function computes all f(sigma) values for all possible vertex subsets, sigma, with size at most max_size.

    Args:
        D (networkx.DiGraph): directed graph
        dictionary (dictionary): dictionary where shortest paths are stored successively
        max_dim (int): maximum size for subsets; default is 3

    Returns:
        all_subsets (list): list of all subsets (list) up to given size
        all_fvals (list): list of all corresponding filtration values (float)
    """

    nodes = sorted(list(D.nodes))
    from itertools import chain, combinations
    all_subsets = list(chain.from_iterable(combinations(nodes, r) for r in range(1,max_size+1)))

    all_fvals = [0 for _ in nodes]
    for sigma in all_subsets[len(nodes):]:
        path = find_f_path(D, sigma, dictionary)
        if path == None:
            all_fvals.append(np.inf)
        else:
            #print(sigma, path)
            #all_fvals.append(len(path)-1)
            all_fvals.append( nx.path_weight(D, path, weight='weight') )
    return all_subsets, all_fvals


def newfiltration_persistence(D, max_dim=1):
    '''
    This function computes persistence via new filtration from a digraph D using Dionysus.

    Args:
        D (networkx.DiGraph): directed graph
        max_dim (int): maximum dimension to compute persistence; default is 1

    Returns:
        dgms (list): list of persistence diagrams
    '''

    shortest_paths_asym = pairs_dictionary(D)
    ## to compute n-dimensional homology we need, (n+1)-dimensional simplices, so subsets of size n+2
    all_subsets, all_fvals = all_fvalues(D, shortest_paths_asym, max_size=max_dim+2)

    f = dio.Filtration()
    for simp, time in zip(all_subsets,all_fvals):
        f.append(dio.Simplex(list(simp), time))
    f.sort()

    dgms = dio.init_diagrams(dio.homology_persistence(f), f)
    return dgms


def plot_dgms(dgms, title=None, filename=None, report_repeats=True, max_dim=1, ax=None, max_val=1, get_dgms=False):
    '''
    Plots persistence diagrams as
    obtained from dionysus.
    '''

    colors = ['red', 'blue', 'green', 'orange', 'brown']

    #if max_dim == None: 
    max_dim = min(max_dim,len(dgms)-1,4)

    births, deaths, dims = [], [], []
    for i,dgm in enumerate(dgms[:max_dim+1]):
        births = births+[pt.birth for pt in dgm]
        deaths = deaths+[pt.death for pt in dgm]
        dims = dims+[i for _ in dgm]
    births = np.array(births)
    deaths = np.array(deaths)

    if report_repeats:
        pts_list = list(zip(births,deaths,dims))

    max_death = max(deaths[deaths<np.inf])
    max_birth = max(births[births<np.inf])
    inf_val = 1.1*max(max_death, max_birth, max_val)

    deaths[ deaths==np.inf ] = inf_val

    if ax == None: ax = plt.gca()

    #fig = plt.figure(figsize=(5.5,5))

    if title != None: ax.title(title)

    # Diagonal
    ax.plot([0,100],[0,100], 'k--', alpha=0.4)
    #plt.plot([0,100],[0,100], 'k--', alpha=0.4)

    ax.set_xticks(range(int(inf_val)+1))
    ax.set_yticks(range(int(inf_val)+1))
    #ax.yticks([inf_val,], minor=True, labels=['\u221e'])

    ax.set_xlim((-0.05*inf_val,inf_val*0.95))
    #ax.xlim((-0.05*max_birth,max_birth*0.95))
    ax.set_ylim((0,inf_val*1.05))
    # Infinity line
    ax.plot([-0.1*inf_val,inf_val*1.5],[inf_val,inf_val], 'k-', alpha=0.3, label='\u221e')

    ax.grid(alpha = 0.2)

    visible_dims = [False]*5
    for i in range(len(dims)):
        ax.plot([births[i],], [deaths[i],], 'o', c=colors[dims[i]], alpha=0.5 )
        visible_dims[dims[i]] = True

    for i,d in enumerate(visible_dims):
        if d: ax.plot([-1,],[-1,], 'o', c=colors[i], label='Dim '+str(i))

    ax.legend(loc='lower right')

    if filename != None: ax.savefig(filename+'.png', dpi=100)

    if report_repeats:
        #pts_list = list(zip(births,deaths))
        count_dict = {i:pts_list.count(i) for i in pts_list}
        counts = np.array(list(count_dict.values()))
        repeats = np.array(list(count_dict.keys()))[counts > 1]
        print('---')
        for r in repeats:
            print('Point', tuple(r), 'shows up', count_dict[tuple(r)], 'times')

    if get_dgms:
        return ax, np.array([births,deaths,dims])

    return ax


def plot_graph_dgms(D, dgms, title=None, report_repeats=True, max_dim=1, max_val=1, spring_iterations=0, with_labels=False, node_size=100, pos=None):
    #fig = plt.figure(figsize=(10,5))
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(6,3))
    if title!=None: plt.suptitle( title )

    #plt.subplot(1,2,1)
    #pos = nx.spring_layout(D,seed=1)
    if spring_iterations>0:
        pos = nx.spring_layout(D, seed=1, k=1, iterations = spring_iterations*len(D))
    elif pos == None:
        n_vert = len(D)
        pos = {i:[np.cos(2*np.pi*i/n_vert),np.sin(2*np.pi*i/n_vert)] for i in D.nodes}
    nx.draw(D, pos=pos, labels={ list(D)[i] : i+1 for i in range(D.number_of_nodes()) }, 
            node_size=node_size, with_labels=with_labels, ax=ax1)

    #plt.subplot(1,2,2)
    plot_dgms(dgms, report_repeats=report_repeats, max_dim=max_dim, max_val=max_val, ax=ax2)

    return fig