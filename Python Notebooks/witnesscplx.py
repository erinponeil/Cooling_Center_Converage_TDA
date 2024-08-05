import gudhi
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import matplotlib
from matplotlib import colors
from sklearn.metrics import pairwise_distances
from geopy import distance
import networkx as nx


def _pairwise_geodesic_distances(landmarks, witnesses):

    D = np.zeros([len(landmarks), len(witnesses)])
    for i in range(len(landmarks)):
        for j in range(len(witnesses)):
            D[i,j] = distance.distance(landmarks[i], witnesses[j]).km

    return D
    

def WitnessComplexFiltration(L,W, verts_immed=True, metric='geodesic'):
    '''
    Computes the filtered simplicial complex via the Witness Complex.
    
    Inputs
    ------
        L: landmark points (these become the vertices of the complex) (np.array)
        W: witness  (np.array)
        verts_immed: if True, vertices are all born at zero, if False, vertices are born at their distance from a witness

    Returns
    -------
        gudhi.SimplexTree() representing filtered simplicial complex (maxdim is 2)
        filtration array saying when vertices merge with other components (i.e. death times)

    '''

    if metric in ['euclidean', 'l1']:
        D = pairwise_distances(L, W, metric='l1')

    else: 
        if metric != 'geodesic':
            print('invalid distance metric. options are euclidean, l1, and geodesic.')
            print('defaulting to geodesic.')

        D = _pairwise_geodesic_distances(L, W)


    filt = np.zeros([len(L), len(L)])
    
    for i in range(len(L)): 

        # vertices all born at 0
        if verts_immed:
            filt[i,i] = 0

        # vertices are born at their distance from a witness
        else:
            filtval = min(D[i,:])
            filt[i,i] = filtval

        # closest witness to vertex i
        am_i = np.argmin(D[i,:])
    
        for j in range(len(L)):
            if j > i:
                # closest witness to vertex j
                am_j = np.argmin(D[j,:])

                # if it's the same witness, the edge is added at the max of dist from j to am_i and i to am_j
                if am_i == am_j:
                    filt[i,j] = max(D[j,am_i], D[i,am_j])

                # if they are different witnesses, its the max of the distance from i to it's closest witness,
                # distance from j to it's closest witness, and 
                else:
                    filt[i,j] = np.amax( [D[i,am_i], D[j,am_j], min(D[j,am_i], D[i,am_j])])

                # symmetric
                filt[j,i] = filt[i,j]
    
    # make simplex tree and create filtration
    st = gudhi.SimplexTree()
    st = st.create_from_array(filt)

    # if it's a valid filtration, add in all triangles to make a clique complex 
    # and return simplex tree
    if not st.make_filtration_non_decreasing():
        st.expansion(2)
        return st, filt

    # If it's not a valid filtration, print a message and return empty simplex tree
    else:
        print('not a filtration')
        return gudhi.SimplexTree(), filt




def GetSublevelSet(st, L, f):
    '''
    Gets list of vertices and edges in sublevelset 

    Inputs
    ------
        st: simplex tree representing filtered witness complex 
        L: landmark points (these become the vertices of the complex) (np.array)
        f: filtration value for sublevelset 

    Returns
    -------
        verts: vertices in sublevelset
        edges: edges in sublevelset
        
    '''

    
    verts = []
    edges = []
    for i in range(len(L)):
        # check if i has filtration value less than f
        if st.filtration([i]) <= f:
            verts.append(L[i,:]) ## add to verts in sublevelset
            
            for j in range(len(L)):
                if j >= i:
                    # check if j has filtration value less than f
                    if st.filtration([j]) <= f:
                        verts.append(L[j,:]) ## add to verts in sublevelset

                        # if i and j are in sublevelset, see if edge between is in sublevelset
                        if st.filtration([i,j]) <= f:
                            edges.append([ L[i,:], L[j,:]]) ## add to edges in sublevelset

    verts = np.array(verts)

    return verts, edges


def GetComponents(filt, seed=24):
    '''
    Gets list of blocks in a component right before it dies

    Inputs
    ------
        filt: filtration array saying when vertices merge with other components (this is an output of WitnessComplexFiltration)

    Returns
    -------
        death: list of death times
        block_list: list of lists of blocks corresponding to each death time 

    e.g. the values in block_list[i] correspond to the indices in CT of the census
    blocks that are in the component that dies at death[i]
        
    '''

    # We need to give our vertices random birth times so that they are distinguishable
    np.random.seed(seed)
    random_births = np.random.rand(np.shape(filt)[0])/100
    rb_dictionary = {i: random_births[i] for i in range(np.shape(filt)[0])}

    # Make a networkx graph with edge weights according to the filtration value
    # The name of the weight is 'score' for no reason other than I'm too lazy to change Gill's code
    G = nx.from_numpy_array(filt, edge_attr='score')
    nx.set_node_attributes(G, rb_dictionary, 'score')
    birth, death, generators, block_list = union_find(G)

    # The only way picking random birth times can impact the result is if we accidentally
    # make it so that edges start being added before all vertices have been added
    if max(birth) > min(death):
        print('Error! Make birth times smaller.')
        return [], []

    return death, block_list



####### The following code was written by Gillian Grindstaff and used here with permission 
def _find(elt,parent):
  while (elt != parent[elt]):
    elt = parent[elt]
  return elt

def _get_children(elt,parent):
  indx = []
  for i,mama in enumerate(parent):
    if mama==elt:
      indx.append(i)
  return indx

def _find_desc(elt,parent):
  desc = [elt]
  children = _get_children(elt,parent)
  for child in children:
    grandbabies = _find_desc(child,parent)
    desc = desc+grandbabies
  return desc


def _union(G,A,B,parent,birth,death,generators,children,block_list,t):
  rootA = _find(A,parent)
  rootB = _find(B,parent)
  if(rootA==rootB):
    return birth, death
  else:
    #elder rule: if A is born first, B dies.
    #We add the birth and death times for B here.
    if G.nodes[rootA]['score']< G.nodes[rootB]['score']:
      birth.append(G.nodes[rootB]['score'])
      parent[rootB] = rootA
      children[rootA] = children[rootA]+children[rootB]
      death.append(t)
      generators.append(rootB)
      block_list.append(children[rootB])
    #if B is born first, A dies
    else:
      parent[rootA] = rootB
      children[rootB] = children[rootB] + children[rootA]
      birth.append(G.nodes[rootA]['score'])
      death.append(t)
      generators.append(rootA)
      block_list.append(children[rootA])
  return _union(G,A,B,parent,birth,death,generators,children,block_list,t)

def union_find(G):
  parent = list(G.nodes)
  birth = []
  death = []
  generators = []
  block_list = []
  children = [[i] for i in G.nodes]
  sorted_by_score = sorted(G.edges(data=True), key=lambda edge_data: edge_data[2]["score"])
  for edge in sorted_by_score:
    t = edge[2]['score']
    v = edge[0]
    w = edge[1]
#    (v,w) = edge.nodes()
    birth,death = _union(G,v,w,parent,birth,death,generators,children,block_list,t)

  bad = []
  #this can be rewritten with np.where more efficiently I think
  for i in range(len(birth)):
    if death[i]-birth[i]==0:
      bad.append(i)

  #remove features with 0 persistence, maybe can be done implicitly
  death = [death[i] for i in range(len(death)) if i not in bad]
  birth = [birth[i] for i in range(len(birth)) if i not in bad]
  generators = [generators[i] for i in range(len(generators)) if i not in bad]
  block_list = [block_list[i] for i in range(len(block_list)) if i not in bad]
  #block_list = [_find_desc(i,parent) for i in tqdm(generators)]
  return birth,death,generators,block_list