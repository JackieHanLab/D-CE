# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 13:52:55 2021

@author: lenovo
"""
import numpy as np
import anndata as ad
import scanpy as sc
import argparse
import numpy as np
import os
#import novosparc
from scipy.spatial.distance import cdist
import numpy as np
# from sklearn import manifold, datasets # not used, remove in next version
from sklearn.neighbors import kneighbors_graph
from scipy.spatial.distance import cdist
from scipy.cluster import hierarchy
from scipy.stats import pearsonr
import ot
import time
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
import os
#import numpy as np
from ot.bregman import sinkhorn
from ot.utils import dist



def tensor_square_loss_adjusted(C1, C2, T):
    """
    Returns the value of \mathcal{L}(C1,C2) \otimes T with the square loss
    function as the loss function of Gromow-Wasserstein discrepancy.

    Where :
        C1 : Metric cost matrix in the source space
        C2 : Metric cost matrix in the target space
        T : A coupling between those two spaces

    The square-loss function L(a,b)=(1/2)*|a-b|^2 is read as :
        L(a,b) = f1(a)+f2(b)-h1(a)*h2(b) with :
            f1(a)=(a^2)/2
            f2(b)=(b^2)/2
            h1(a)=a
            h2(b)=b

    Parameters
    ----------
    C1 : ndarray, shape (ns, ns)
         Metric cost matrix in the source space
    C2 : ndarray, shape (nt, nt)
         Metric costfr matrix in the target space
    T : ndarray, shape (ns, nt)
         Coupling between source and target spaces

    Returns
    -------
    tens : ndarray, shape (ns, nt)
           \mathcal{L}(C1,C2) \otimes T tensor-matrix multiplication result
    """

    C1 = np.asarray(C1, dtype=np.float64)
    C2 = np.asarray(C2, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)

    def f1(a):
        return (a**2) / 2

    def f2(b):
        return (b**2) / 2

    def h1(a):
        return a

    def h2(b):
        return b

    tens = -np.dot(h1(C1), T).dot(h2(C2).T) 
    tens -= tens.min()

    return tens



def gromov_wasserstein_adjusted_norm(cost_mat, C1, C2, alpha_linear,p, q, loss_fun, epsilon,
                                     max_iter=1000, tol=1e-9, verbose=False, log=False):
    """
    Returns the gromov-wasserstein coupling between the two measured similarity matrices

    (C1,p) and (C2,q)

    The function solves the following optimization problem:

    .. math::
        \GW = arg\min_T \sum_{i,j,k,l} L(C1_{i,k},C2_{j,l})*T_{i,j}*T_{k,l}-\epsilon(H(T))

        s.t. \GW 1 = p

             \GW^T 1= q

             \GW\geq 0

    Where :
        M  : cost matrix in sourceXtarget space
        C1 : Metric cost matrix in the source space
        C2 : Metric cost matrix in the target space
        p  : distribution in the source space
        q  : distribution in the target space
        L  : loss function to account for the misfit between the similarity matrices
        H  : entropy

    Parameters
    ----------
    M : ndarray, shape (ns, nt)
         Cost matrix in the sourceXtarget space
    C1 : ndarray, shape (ns, ns)
         Metric cost matrix in the source space
    C2 : ndarray, shape (nt, nt)
         Metric costfr matrix in the target space
    p :  ndarray, shape (ns,)
         distribution in the source space
    q :  ndarray, shape (nt,)
         distribution in the target space
    loss_fun :  string
        loss function used for the solver either 'square_loss' or 'kl_loss'
    epsilon : float
        Regularization term >0
    max_iter : int, optional
       Max number of iterations
    tol : float, optional
        Stop threshold on error (>0)
    verbose : bool, optional
        Print information along iterations
    log : bool, optional
        record log if True

    Returns
    -------
    T : ndarray, shape (ns, nt)
        coupling between the two spaces that minimizes :
            \sum_{i,j,k,l} L(C1_{i,k},C2_{j,l})*T_{i,j}*T_{k,l}-\epsilon(H(T))
    """

    C1 = np.asarray(C1, dtype=np.float64)
    C2 = np.asarray(C2, dtype=np.float64)
    cost_mat = np.asarray(cost_mat, dtype=np.float64)

    T = np.outer(p, q)  # Initialization

    cpt = 0
    err = 1
    
    try:
            cost_mat_norm = cost_mat/ cost_mat.max()
    except:
            cost_mat_norm = cost_mat
            
    if alpha_linear == 1:
        T = sinkhorn(p, q, cost_mat_norm, epsilon)
    else:        
        while (err > tol and cpt < max_iter):
            
            Tprev = T

            if loss_fun == 'square_loss':
                tens = tensor_square_loss_adjusted(C1, C2, T)

            tens_all = (1-alpha_linear)*tens + alpha_linear*cost_mat_norm
            T = sinkhorn(p, q, tens_all, epsilon)
        
            if cpt % 10 == 0:
            # We can speed up the process by checking for the error only all
            # the 10th iterations
                err = np.linalg.norm(T - Tprev)

                if log:
                    log['err'].append(err)

                if verbose:
                    if cpt % 200 == 0:
                        print('{:5s}|{:12s}'.format(
                                'It.', 'Err') + '\n' + '-' * 19)
                        print('{:5d}|{:8e}|'.format(cpt, err))

            cpt += 1

    if log:
        return T, log
    else:
        return T



def load_target_space(path, cells_selected=None, is_2D=True):
    locations = np.loadtxt(path, skiprows=1)
    if is_2D:
        locations = locations[:, [0, 2]]
    if cells_selected is not None:
        locations = locations[cells_selected, :]
    return locations


def load_data(path, dtype='dge'):
    if dtype == 'dge':
        dataset = ad.read_text(path)

    elif dtype == '10x':
        dataset = sc.read_10x_mtx(path,  var_names='gene_symbols',  cache=True)

    return dataset

def reconstruct_structure(dataset, locations, insitu_matrix=None,exp_matrix = None):
    dataset=dataset
    dge = dataset.X
    num_cells = len(dataset.obs)
    num_locations = locations.shape[0]
    locations = locations
    insitu_matrix=insitu_matrix
    exp_matrix = exp_matrix
    if exp_matrix is None:
        cost_marker_genes = np.ones((num_cells, num_locations))
    else:
        print('marker used')
        cost_marker_genes = cdist(exp_matrix/np.amax(exp_matrix),insitu_matrix/np.amax(insitu_matrix))
    num_neighbors_target= 7
    A_locations = kneighbors_graph(locations, num_neighbors_target, mode='connectivity', include_self=True)
    sp_locations = dijkstra(csgraph = csr_matrix(A_locations), directed = False,return_predecessors = False)
    sp_locations_max = np.nanmax(sp_locations[sp_locations != np.inf])
    sp_locations[sp_locations > sp_locations_max] = sp_locations_max #set threshold for shortest paths
    num_neighbors_source = 7 # number of neighbors for nearest neighbors graph at source
    A_expression = kneighbors_graph(dge, num_neighbors_source, mode='connectivity', include_self=True)
    sp_expression = dijkstra(csgraph = csr_matrix(A_expression), directed = False, return_predecessors = False) 
    sp_expression_max = np.nanmax(sp_expression[sp_expression != np.inf])
    sp_expression[sp_expression > sp_expression_max] = sp_expression_max #set threshold for shortest paths

    # Set normalized cost matrices based on shortest paths matrices at target and source spaces
    cost_locations = sp_locations / sp_locations.max()
    cost_locations -= np.mean(cost_locations)
    cost_expression = sp_expression / sp_expression.max()
    cost_expression -= np.mean(cost_expression)
    costs = {'expression':cost_expression,'locations': cost_locations,'markers': cost_marker_genes}
    p_locations = ot.unif(num_locations)
    p_expression = ot.unif(num_cells)
    epsilon=5e-4
    if exp_matrix is None:
        gw = gromov_wasserstein_adjusted_norm(cost_marker_genes, cost_expression, cost_locations,
												  0, p_expression, p_locations,
												  'square_loss', epsilon=epsilon, verbose=True)
    else:
        gw = gromov_wasserstein_adjusted_norm(cost_marker_genes, cost_expression, cost_locations,
												  0.5, p_expression, p_locations,
												  'square_loss', epsilon=epsilon, verbose=True)
    coordnew = np.dot(gw,locations)
    return gw,coordnew
    
def mapping(gw,target):
    gw1=gw
    nodeindex=list(range(0,gw.shape[1]))
    cellindex=list(range(0,gw.shape[0]))
    nodeindex2=nodeindex
    for a in range(0,gw.shape[1]):
        matchage2=gw1[cellindex,nodeindex2]
        maxindex=np.where(gw1==np.max(matchage2))
        for b in range(0,len(maxindex[0])):
            if maxindex[0][b] in cellindex and maxindex[1][b] in nodeindex2:
                index=[maxindex[0][b],maxindex[1][b]]
                break
        weight=sum(gw1[:,index[1]])
        gw1[:,index[1]]=0
        gw1[index[0],index[1]]=weight
        nodeindex2.remove(maxindex[1][b])
        cellindex.remove(maxindex[0][b])
    for a in range(0,gw.shape[1]):
        gw1[a,:]=gw1[a,:]/sum(gw1[a,:])
    coordnew=np.dot(gw1,target)
    return coordnew


parser = argparse.ArgumentParser()
parser.add_argument('-c','--coord_path',help='reconstructed coordinates using D-CE',dest='coord')
parser.add_argument('-r','--ref_path',help='expression data of reference',dest='ref_matrix',default=None)
parser.add_argument('-e','--exp_path',help='expression matrix',dest='exp_matrix',default=None)
parser.add_argument('-t','--target_space_path',help='coordinate of target space',dest='target_matrix',default=None)
args = parser.parse_args()

dataset = np.loadtxt(args.coord)
locations = load_target_space(args.target_matrix, is_2D=False)

exp_path = args.exp_matrix
if not exp_path:
    adata = ad.AnnData(dataset)
    result=reconstruct_structure(adata,locations)
else:
    adata = ad.AnnData(dataset)
    insitu_matrix = np.loadtxt(args.ref_matrix)
    exp = np.loadtxt(args.exp_matrix)
    result = reconstruct_structure(adata,locations,insitu_matrix,exp)
coordnew=mapping(result[0],locations)
np.savetxt('coordnew.txt',coordnew,fmt='%.4e')
