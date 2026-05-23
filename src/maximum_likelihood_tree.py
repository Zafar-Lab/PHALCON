# IMPORTS ###########################################################################################################
import sys
import random
import numpy as np
import pandas as pd
from tqdm import tqdm
import time
import copy
from ete3 import Tree
from src import NNI
from src import SPR
#####################################################################################################################




def getAltAlleleCounts(arr):
# Converts read count strings (e.g. "0,1,2,1") to numpy array (e.g. [0,1,2,1])
    return np.array(list(map(int,arr.strip().split(","))))




def chunk_ML(splits, inverted_splits, ones, zeros,split_names):
# Calculate the maximum likelihood for each site

    tmp = np.matmul(splits,ones)+np.matmul(inverted_splits,zeros)
    pd.DataFrame(tmp).to_csv("tmp_matrix.tsv",sep='\t',index=False,header=False)
    pd.DataFrame(split_names).to_csv("genotype_name_tmp.tsv",sep='\t',index=False,header=False)
    par_case = 0
    for i in range(split_names.shape[0]):
       if split_names[i][-1] == 'Par':
          par_case+=1

    loss_case = 0
    for i in range(split_names.shape[0]):
       if split_names[i][-1] == 'Loss':
          loss_case+=1


    priors = np.zeros((splits.shape[0],ones.shape[1]))
    for i in range(split_names.shape[0]):
       if split_names[i][-1] == 'Par':
          priors[i] = np.log(  np.exp(-700) * 0.001 * (1/par_case)  ) 
       elif split_names[i][-1] == 'Loss':
          priors[i] = np.log(  np.exp(-700) * 0.001 * (1/loss_case) ) 

    tmp = tmp + priors
    split_indexes = np.argmax(tmp,axis=0)
    genotypes = splits[split_indexes]
    split_names = split_names[split_indexes]
    L = np.sum(np.max(tmp, axis=0))
    return (L, genotypes, split_names)




def heterozygous(ccm, tree, rc, names, split_list, split_counts, split_names,clusters):
   for node in tree.traverse("preorder"):
        if node.is_leaf():
          mutated_names = set([node.name])
        else:
          mutated_names = set(node.get_leaf_names(is_leaf_fn=None))
        tmp_arr = np.zeros(rc.shape[1])
        for name in names:
           if name in list(mutated_names):
              for cell in ccm[name]:
                tmp_arr[cell] = 1
        tmp_arr = tmp_arr.tolist()
        tmp_copy = tmp_arr.copy()
        if tmp_arr not in split_list:
            split_list.append(tmp_arr)
            split_counts.append(tmp_arr.count(1))
            if tmp_arr.count(1) == len(tmp_arr):
               tmp_copy.append('Clonal')
            else:
               tmp_copy.append('Het')
            split_names.append(tmp_copy)
   return split_counts, split_list, split_names




def parallel(ccm, tree, rc, names, split_list, split_counts,split_names,clusters):
   m = []  #list to hold mutated names (includes duplicates)
   for node1 in tree.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  # descendants of current node
     ans = list( h for h in node1.get_ancestors() )  # anscestors of current node
     for node2 in tree.traverse("preorder"):
      if len(ans)==0:
       break
      else:
       if node1 != node2 and node2 not in ans:
         if node2 not in des and node1.up!=node2.up:
           if node1.is_leaf():
             mutated_names1 = set([node1.name])
           else:
             mutated_names1 = set(node1.get_leaf_names(is_leaf_fn=None))
           if node2.is_leaf():
             mutated_names2 = set([node2.name])
           else:
             mutated_names2 = set(node2.get_leaf_names(is_leaf_fn=None))
           m.append( mutated_names1.union( mutated_names2))
   mutated_names = []
   for i in m:
     if i not in mutated_names:
       mutated_names.append(i)    # List of all possibilities of parallel mutation (no duplicate cases)
   for j in mutated_names:
     not_mutated_names = clusters - j
     tmp_arr = np.zeros(rc.shape[1])
     for name in names:
       if name in list(j):
         for cell in ccm[name]:
           tmp_arr[cell] = 1
     tmp_arr = tmp_arr.tolist()
     tmp_copy = tmp_arr.copy()
     if tmp_arr not in split_list:
       split_list.append(tmp_arr)
       split_counts.append(tmp_arr.count(1))
       tmp_copy.append('Par')
       split_names.append(tmp_copy)
   return split_counts, split_list, split_names



def loss(ccm, tree, rc, names, split_list, split_counts,split_names, clusters):
   m = []  #another list to hold mutated names (includes duplicates)
   for node1 in tree.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  #list of descendants of current node
     ans = list( h for h in node1.get_ancestors() )  #list of anscestors of current node
     for node2 in tree.traverse("preorder"):
      if len(des)==0:
       break
      else:
       if node1 != node2 and node2 not in ans:
         if node2 in des:
           if node1.is_leaf():
             mutated_names1 = set([node1.name])
           else:
             mutated_names1 = set(node1.get_leaf_names(is_leaf_fn=None))
           if node2.is_leaf():
             non_mutated_names2 = set([node2.name])
           else:
             non_mutated_names2 = set(node2.get_leaf_names(is_leaf_fn=None))
           m.append( mutated_names1 - non_mutated_names2)
   mutated_names = []
   for i in m:
    if i not in mutated_names:
     mutated_names.append(i)
   for j in mutated_names:
     not_mutated_names = clusters - j
     tmp_arr = np.zeros(rc.shape[1])
     for name in names:
       if name in list(j):
         for cell in ccm[name]:
           tmp_arr[cell] = 1
     tmp_arr = tmp_arr.tolist()
     tmp_copy = tmp_arr.copy()
     if tmp_arr not in split_list:
       split_list.append(tmp_arr)
       split_counts.append(tmp_arr.count(1))
       tmp_copy.append('Loss')
       split_names.append(tmp_copy)
   return split_counts, split_list, split_names
   


def ML(ccm, tree, rc, names, one_Ls, zero_Ls):  
  # ccm : cluster cell map , 
  # rc : Read count matrix (site x cell), 
  # one_Ls : likelihood matrix of mutation occuring (sites x cells)
  # zero_Ls : likelihood of mutation not occuring (sites x cells) 
  clusters = copy.copy(names)
  clusters = set(clusters)
  split_list = []
  split_counts = []
  split_names = []

  split_counts, split_list, split_names = heterozygous(ccm, tree, rc, names, split_list, split_counts,split_names, clusters)
  split_counts, split_list, split_names = parallel(ccm, tree, rc, names, split_list, split_counts,split_names, clusters)
  split_counts, split_list, split_names = loss(ccm, tree, rc, names, split_list, split_counts,split_names, clusters)

  tmp_arr = np.zeros(rc.shape[1])  # The case for no mutation at all
  tmp_arr = tmp_arr.tolist()
  tmp_copy = tmp_arr.copy()
  split_list.append(tmp_arr) 
  split_counts.append(tmp_arr.count(1))
  tmp_copy.append('None')
  split_names.append(tmp_copy)

  split_list = np.array(split_list)
  split_names = np.array(split_names,dtype='object')

  indexes = np.argsort(split_counts)
  splits_sorted = split_list[indexes,:]
  split_names = split_names[indexes,:]
  new_genotypes = np.zeros((rc.shape[0], rc.shape[1]))
  L_wg = np.float128(0)
  s = time.time()
  (L_wg, new_genotypes,split_names) = chunk_ML(splits=splits_sorted, inverted_splits=1-splits_sorted, ones=np.float16(one_Ls), zeros=np.float16(zero_Ls),split_names = split_names)
  new_genotypes = new_genotypes.T
  return (new_genotypes,L_wg,split_names)

  

  

def ML_initialization(lc, one_Ls, zero_Ls):
# Initialisation of likelihood matrices (no mutation likelihood matrix and mutation likelihood matrix)
# Initialisation of the likelihood value
    init_L = np.float128(0)
    num = lc.shape[1]
    pos = lc.shape[0]
    lc = lc.replace(0,0.00000001)
    lc = lc.replace(1,0.99999999)
    for i in tqdm(range(pos),desc="Initializing Log likelihoods..."):
        for j in range(num):
            one_Ls[i][j] = np.float16(np.log(lc.iloc[i][j]))
            zero_Ls[i][j] = np.float16(np.log(1-lc.iloc[i][j]))
    one_Ls = one_Ls.T
    zero_Ls = zero_Ls.T
    return (init_L, one_Ls, zero_Ls)



def findBestLikelihood(iterations,n_cluster,data_df,one_Ls,zero_Ls,final_labels,lklhd_computed): 
# Find best likelihoods among a set of chosen number of clusters

    names=[]
    for i in range(n_cluster):
        name = "c"+str(i+1)
        names.append(name)
    cluster_cell_map={}
    for name in names:
        ind = int(name[1:])
        cluster_cell_map[name] = [i for i in range(len(final_labels)) if final_labels[i] == ind-1]
  

    mut_tree = Tree()
    mut_tree.populate(len(names),names)

    ete_nj_init = Tree()
    ete_nj_init.populate(len(names),names)

    best_L = float("-inf")
    (new_mat, Likelihood, split_names) = ML(ccm=cluster_cell_map,tree=ete_nj_init, rc=data_df.T, names=names, one_Ls=one_Ls, zero_Ls=zero_Ls)

    if Likelihood>best_L:
        best_L = Likelihood

    n_iterations = iterations
    best_Ls = [best_L]
    stack  = [ete_nj_init]
    top_ids = set()
    top_ids.add(ete_nj_init.get_topology_id())
    best_mat = new_mat
    best_tree = ete_nj_init
    best_geno_name = split_names
    best_trees = [best_tree]
    best_geno_names = [split_names]
    for it in range(n_iterations):
        print("Iteration # ",it)
        Ls = []
        ts = []
        mats = []
        genos = []
        for item_ in stack:
            tree_list = NNI.Main(in_tree=item_, N=n_cluster)
            for tree_ in tree_list:
                if tree_.get_topology_id() not in top_ids:
                    top_ids.add(tree_.get_topology_id())
                    (mat_, Likelihood, split_names) = ML(ccm=cluster_cell_map,tree=tree_, rc=lklhd_computed, names=names, one_Ls=one_Ls, zero_Ls=zero_Ls)
                    if Likelihood > best_L:
                        best_L = Likelihood
                        best_mat = mat_
                        best_tree = tree_
                        best_geno_name = split_names
                        Ls=[]
                        mats=[]
                        ts=[]
                        genos = []
                        Ls.append(Likelihood)
                        mats.append(mat_)
                        ts.append(tree_)
                        genos.append(split_names)
                
            tree_list = SPR.Main(in_tree=item_, N=n_cluster, N_dest=n_cluster)
            for tree_ in tree_list:
                if tree_.get_topology_id() not in top_ids:
                    top_ids.add(tree_.get_topology_id())
                    (mat_, Likelihood, split_names) = ML(ccm=cluster_cell_map,tree=tree_, rc=lklhd_computed, names=names, one_Ls=one_Ls, zero_Ls=zero_Ls)
                    if Likelihood > best_L:
                        best_L = Likelihood
                        best_mat = mat_
                        best_tree = tree_
                        best_geno_name = split_names
                        Ls=[]
                        mats=[]
                        ts=[]
                        genos=[]
                        Ls.append(Likelihood)
                        mats.append(mat_)
                        ts.append(tree_)
                        genos.append(split_names)
                
        max_ = float("-inf")
        if len(Ls)!=0:
            stack = ts
        else:
            print("no more better proposed trees")
            print("terminating the search")
            break
    return best_L,best_mat,best_tree, best_geno_name 
                 





