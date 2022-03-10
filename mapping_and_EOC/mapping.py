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
import novosparc
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
    if not exp_matrix:
        cost_marker_genes = np.ones((num_cells, num_locations))
    else:
        cost_marker_genes = cdist(exp_matrix/np.amax(exp_matrix),insitu_matrix/np.amax(insitu_matrix))
    num_neighbors_target= 5
    A_locations = kneighbors_graph(locations, num_neighbors_target, mode='connectivity', include_self=True)
    sp_locations = dijkstra(csgraph = csr_matrix(A_locations), directed = False,return_predecessors = False)
    sp_locations_max = np.nanmax(sp_locations[sp_locations != np.inf])
    sp_locations[sp_locations > sp_locations_max] = sp_locations_max #set threshold for shortest paths
    num_neighbors_source = 5 # number of neighbors for nearest neighbors graph at source
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
    p_locations, p_expression = novosparc.rc.create_space_distributions(num_locations, num_cells)
    epsilon=5e-4
    if not exp_matrix:
        gw = novosparc.rc._GWadjusted.gromov_wasserstein_adjusted_norm(cost_marker_genes, cost_expression, cost_locations,
												  0, p_expression, p_locations,
												  'square_loss', epsilon=epsilon, verbose=True)
    else:
        gw = novosparc.rc._GWadjusted.gromov_wasserstein_adjusted_norm(cost_marker_genes, cost_expression, cost_locations,
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
args=parser.parse_args(['-c', 'coordznew.txt', '-t', 'geometry.txt'])
#ref_matrix = load_data(args.ref_matrix).T
dataset = np.loadtxt(args.coord)
locations = load_target_space(args.target_matrix, is_2D=False)
locations = locations[0:3039,:]
exp_path = args.exp_matrix
if not exp_path:
    adata = ad.AnnData(dataset)
    result=reconstruct_structure(adata,locations)
else:
    insitu_matrix = np.loadtxt(args.ref_matrix)
    exp = np.loadtxt(args.exp_matrix)
    result = reconstruct_structure(adata,target_space_path,insitu_matrix,exp)
coordnew=mapping(result[0],locations)
np.savetxt('coordnew.txt',coordnew,fmt='%.4e')
