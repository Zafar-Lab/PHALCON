import sys
import random
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import betabinom
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import scipy
from scipy.sparse import csgraph
from scipy.sparse.linalg import eigsh,eigs
from scipy.linalg import eigvalsh
from sklearn.manifold import TSNE
import seaborn as sns
from sklearn.impute import KNNImputer
import time
import loompy
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import f1_score, recall_score, precision_score
from sklearn.neighbors import NearestCentroid
import copy
import os
import csv
import math
from decimal import *
plt.switch_backend('agg')
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
from sklearn.metrics.cluster import adjusted_rand_score

def indexToChar(index):
    if (index == 0):
        return 'A'
    elif (index == 1):
        return 'C'
    elif (index == 2):
        return 'G'
    elif (index == 3):
        return 'T'

def charToIndex(char):
    if (char == 'A'):
        return 0
    elif (char == 'C'):
        return 1
    elif (char == 'G'):
        return 2
    elif (char == 'T'):
        return 3


with loompy.connect('/home/priya/Downloads/Final Stage/Mission_Bio_looms/MB_TN1_original.cells.loom') as ds:
        print('ds layers: ',ds.layers)
        print("ds shape: ",ds.shape)
        gt=ds[''][:]
        print("gt shape: ",gt.shape)
        dp=ds.layers.DP[:]      
        print("dp shape: ",dp.shape)
        gq = ds.layers.GQ[:]
        ad=ds['AD'][:]
        
        print("ad shape: ",ad.shape)
        ro=ds['RO'][:]
        print("ro shape: ",ro.shape)
#         gq=ds['GQ'][:]
#         print("gq shape: ",gq.shape)
        rows = ds.shape[0]
        cols = ds.shape[1]
        print("rows: ",rows)
        print("cols: ",cols)
        snp=pd.DataFrame({'CHROM':ds.ra['CHROM'],'POS':ds.ra['POS'],'REF':ds.ra['REF'],'ALT':ds.ra['ALT']})
        gt_df = pd.DataFrame(gt, index=None)
        gt_df.index = "chr"+snp["CHROM"].map(str) + "_"+snp["POS"].map(str)+"_"+snp["REF"].map(str)+"_"+snp["ALT"].map(str)



all_pos_list = list(gt_df.index)
ind_to_be_filtered_out = []
for ind,pos in enumerate(all_pos_list):
    pos_arr = pos.split('_')
    if len(pos_arr[2])>1 or len(pos_arr[3])>1:
        ind_to_be_filtered_out.append(ind)
print(len(ind_to_be_filtered_out))
gq = np.delete(gq,ind_to_be_filtered_out,axis=0)
dp = np.delete(dp,ind_to_be_filtered_out,axis=0)
ad = np.delete(ad,ind_to_be_filtered_out,axis=0)
ro = np.delete(ro,ind_to_be_filtered_out,axis=0)


mpileupFileName = "output_AML_nodups.mpileup"

rcFileName = "readCounts_MBTN1_nodups.tsv"

quality_file = "gq_file_TN1.tsv"

all_cols = gt_df.index.to_numpy()
all_cols = np.delete(all_cols,ind_to_be_filtered_out,axis=0)
print('len(all_cols): ',len(all_cols))

all_rows = gt_df.columns.to_numpy()
print('len(all_rows): ',len(all_rows))

mpileupFile = open(mpileupFileName,'w')

rcFile = open(rcFileName,'w')
gqFile = open(quality_file,'w')
pos_ind=0
c=0
prev = ''
for col in tqdm(all_cols,desc="Constructing mpileup file from loom file.."):
    pos_details = col.split('_')
    chr_name = pos_details[0]
    position = pos_details[1]
    if prev==position:
        c+=1
        continue
    else:
        prev = position
    ref_allele = pos_details[2]
    alt_allele = pos_details[3]
    if alt_allele=='' or alt_allele=='*':
        c+=1
        continue
    mpileupFile.write(chr_name+"\t" + position +"\t" + ref_allele+ "\t")
    rcFile.write(chr_name+"\t" + position +"\t" + ref_allele+ "\t"+alt_allele)
    gqFile.write(position)
    for r in range(len(all_rows)):

        alt_depth = ad[c][r]
        coverage = dp[c][r]
        ref_depth = ro[c][r]
        gq_qual = gq[c][r]
        readString = ""

        if coverage == 0:
            mpileupFile.write("0\t"+"*\t*")
            rcFile.write("\t"+"0,0,0,0")
            gqFile.write("\t"+"0")
        else:
            baseQualityStr = "I" * coverage
            mpileupFile.write(str(coverage)+"\t")
            readString += "." * int(ref_depth)
            readString += alt_allele * int(alt_depth)
            counts = np.zeros(4)
            counts[charToIndex(ref_allele)] = int(ref_depth)
            counts[charToIndex(alt_allele)] = int(alt_depth)
            if coverage > ref_depth+alt_depth:
                diff = coverage-(ref_depth+alt_depth)
                other_allele = indexToChar(random.randint(0,3))
                while other_allele == ref_allele or other_allele == alt_allele: 
                    other_allele = indexToChar(random.randint(0,3))
                readString += other_allele * int(diff)
                counts[charToIndex(other_allele)] = int(diff)
            mpileupFile.write(readString+"\t")
            mpileupFile.write(baseQualityStr)
            rcFile.write('\t')
            rcFile.write(str(int(counts[0]))+','+str(int(counts[1]))+','+str(int(counts[2]))+','+str(int(counts[3])))
            gqFile.write("\t"+str(gq_qual))
        if r!=len(all_rows)-1:
            mpileupFile.write("\t")
    mpileupFile.write("\n")   
    rcFile.write('\n')
    gqFile.write('\n')
    c+=1
    
mpileupFile.close()
rcFile.close()
gqFile.close()

