#!/usr/bin/env python3

# IMPORTS ###########################################################################################################
import sys
import random
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import betabinom
from sklearn.cluster import SpectralClustering
from scipy.sparse import csgraph
from sklearn.neighbors import kneighbors_graph
from sklearn.impute import KNNImputer
import time
import copy
import os
from scipy import linalg
import reload 
import importlib 
sys.path.insert(0, os.path.abspath(".")) 
import filtering, maximum_likelihood_tree, write_vcf
importlib.reload(filtering)
importlib.reload(maximum_likelihood_tree)
importlib.reload(write_vcf)
import scanpy as sc
from sklearn.metrics import f1_score, recall_score, precision_score,silhouette_score
import seaborn as sns
import matplotlib.pyplot as plt

#####################################################################################################################
class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def add_argument(self, action):
        if action.help is not argparse.SUPPRESS:
            parts = []
            if action.option_strings:
                parts.append(', '.join(action.option_strings))  
            if action.help:
                parts.append(action.help)
            self._add_item(self._format_action, [action, ' '.join(parts)])

argParser = argparse.ArgumentParser(
    prog='phalcon',
    description="PHALCON: Phylogeny-aware variant calling from large-scale single-cell panel sequencing datasets",
    formatter_class=CustomHelpFormatter
)


argParser = argparse.ArgumentParser(prog='PROG')
argParser.add_argument('-i','--inReadCountFileName', type=str, help = 'Input read count file',default = '/sample_data/sample_readcount_file.tsv')  # CHANGE HERE
argParser.add_argument('-g', '--inQualityFileName', type=str, help="Genotype quality file",default = '/sample_data/sample_gq_file.tsv')
argParser.add_argument('-o', '--outputPrefix', type=str, help = "Output prefix", default='out_')
argParser.add_argument('-r', '--minReadDepth', type=int, help = "Read depth threshold (Default = 5)",default=5)
argParser.add_argument('-q','--minGenotypeQuality',type=int, help = "Genotype quality threshold (Default = 30)",default=30)
argParser.add_argument('-a', '--minAltAlleleFrequency', type=float, help = "Alternate Frequency Threshold (Default = 0.2)", default=0.2)
argParser.add_argument('-v', '--minVarThreshold', type=float, default=0.5)  
argParser.add_argument('-m', '--minMutantFraction', type=float, default=0.004)
argParser.add_argument('-e', '--eigenValueThreshold', type=float, default=0.5) 
argParser.add_argument('-p', '--maxDecreaseOfClusterCount', type=int, default=0)
argParser.add_argument('-l', '--maxIncreaseOfClusterCount', type=int, default=0)
argParser.add_argument('-c','--clusteringTechnique',type=str,default='spectral',help='Options: "spectral", "leiden"')
argParser.add_argument('-gq', '--enableGenotypeQualityFilter', help = "Put 1 to enable the genotype quality filter (Default = 0)",type=int, default=0)
argParser.add_argument('-b','--treeIterations',type=int,default=30)
argParser.add_argument('-s', '--seed', type=int, default = 12099)
argParser.add_argument("-min","--minDist",type=float,default=0.2)


args = argParser.parse_args()


iterations = args.treeIterations
countFileName = args.inReadCountFileName
qualityFileName = args.inQualityFileName
outputPrefixName = args.outputPrefix
seed=args.seed
eigenValueThreshold = args.eigenValueThreshold
lt = args.maxDecreaseOfClusterCount
gt = args.maxIncreaseOfClusterCount
isQualityFilteringOn = args.enableGenotypeQualityFilter
inferredGenotypeFileName=outputPrefixName+'inferred_genotypes.tsv'
minReadDepth = args.minReadDepth
minAltAlleleFreq = args.minAltAlleleFrequency
varThreshold = args.minVarThreshold
minMutFraction = args.minMutantFraction
qualityValue = args.minGenotypeQuality
clustering = args.clusteringTechnique
min_dist_defined = args.minDist


def indexToChar(index):
# Converts index to character
    if (index == 0):
        return 'A'
    elif (index == 1):
        return 'C'
    elif (index == 2):
        return 'G'
    elif (index == 3):
        return 'T'


def charToIndex(char):
# Converts character to index
    if (char == 'A'):
        return 0
    elif (char == 'C'):
        return 1
    elif (char == 'G'):
        return 2
    elif (char == 'T'):
        return 3



def charToNum(t):
# Returns int value of an integer string
    if t.isdigit():
        return int(t)
    return t



def getAltAlleleCounts(arr):
# Converts read count strings (e.g. "0,1,2,1") to numpy array (e.g. [0,1,2,1])
    return np.array(list(map(int,arr.strip().split(","))))





def generate_graph_laplacian(df, nn):
# Generate Graph Laplacian from data 
    connectivity = kneighbors_graph(X=df, n_neighbors=nn, mode='connectivity')
    adjacency_matrix_s = (1/2)*(connectivity + connectivity.T)
    graph_laplacian_s = csgraph.laplacian(csgraph=adjacency_matrix_s, normed=True)#unnormalized laplacian
    graph_laplacian = graph_laplacian_s.toarray()
    return graph_laplacian



def compute_spectrum_graph_laplacian(graph_laplacian):
# Compute eigenvalues and eigenvectors and project them onto the real numbers
    eigenvals, eigenvcts = linalg.eig(graph_laplacian)
    eigenvals = np.real(eigenvals)
    eigenvcts = np.real(eigenvcts)
    return eigenvals, eigenvcts



random.seed(seed)
np.random.seed(seed)


getAltAlleleCounts =np.frompyfunc(getAltAlleleCounts, 1, 1)  


# # # # # # # # # # # # # # # # # # # # # # FINAL ALGORITHM # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
start_time = time.time()
df = pd.read_csv(countFileName, sep="\t",header=None)
df = df.to_numpy()
print("Dimensions of the input matrix: ",df.shape)


# # # # # # # # # # # # # # # # # # # # # # FILTERS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
filters_start_time = time.time()
muts = df.shape[0]
cells = df.shape[1]-4
alt_alleles = []
retain_indices = []

for i in tqdm(range(muts),desc='Computing Alternate alleles...'):
    if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
        allele_counts = np.zeros((cells,4),dtype=int)
        ref = df[i][2]
        ref_ind = charToIndex(ref)
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        allele_counts = np.sum(counts,axis=0)
        ref_count = allele_counts[ref_ind]
        allele_counts[ref_ind] = -1
        if np.max(allele_counts)<=0:
            continue
        retain_indices.append(i)
        alt_ind = np.argmax(allele_counts)
        alt = indexToChar(alt_ind)
        alt_alleles.append(alt)
    elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
        allele_counts = np.zeros((cells,2),dtype=int)
        ref = df[i][2]
        alt = df[i][3]
        ref_ind = 0
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        allele_counts = np.sum(counts,axis=0)
        ref_count = allele_counts[0]  # reference is present at index 0 in case of indels
        allele_counts[ref_ind] = -1
        if np.max(allele_counts)<=0:
            continue
        retain_indices.append(i)
        alt_alleles.append(alt)
    else:
       print("Problem occured while computing alternate alleles at pos ",i )


df = df[retain_indices,:]
df_copy = copy.copy(df)
alt_alleles_copy = copy.copy(alt_alleles)



# .................................Quality Filter..............................................................................................
if isQualityFilteringOn:
    df_gq = pd.read_csv(qualityFileName,sep='\t',header=None)
    gq = df_gq.drop(columns=df_gq.columns[0],axis=1).to_numpy()
    gq = gq[retain_indices,:]
    filtering.genoQualityFilter(df,gq,qualityValue)

# .................................Read depth Filter..............................................................................................
filtering.readDepthFilter(df,minReadDepth)

# .................................Calculating likelihood matrix.............................................................................
muts = df.shape[0]
cells = df.shape[1]-4
lklhd_df = pd.DataFrame(index=range(muts), columns = range(cells))
alt_pos = 0
lklhd_df = lklhd_df.to_numpy()

print('Computing likelihoods... ')
for i in tqdm(range(muts),desc='Computing likelihoods... '):
    alt = alt_alleles[alt_pos]
    alt_pos += 1
    counts = getAltAlleleCounts(df[i][4:])
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)

    
    if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
       altIndex = charToIndex(alt)
       alt_depth = counts[:,altIndex]
    elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
       alt_depth = counts[:,1]  # because in case of indels, the alternate depth is always present at index 1 (ref:alt)

    alpha = (-0.000027183*coverage)+0.068567471
    beta = (0.007454388*coverage)+2.367486659
    w = (0.000548761*coverage)+0.540396786
    alpha_1 = (0.057378844*coverage)+0.669733191
    alpha_2 = (0.003233912*coverage)+0.399261625

    l0 = betabinom.pmf(alt_depth, coverage, alpha, beta)
    l1 = (w*betabinom.pmf(alt_depth, coverage, alpha_1, alpha_1)) + ((1-w)*betabinom.pmf(alt_depth, coverage, alpha_2, alpha_2))
    l2 = betabinom.pmf(alt_depth, coverage, beta, alpha)
    
    l1 = np.where(l1<0,betabinom.pmf(alt_depth, coverage, alpha_1, alpha_1),l1)

    l0 = l0/(l0+l1+l2)
    lklhd_df[i] = np.where(coverage==0,0,1-l0)

print("-------After computing likelihood matrix-------")


# .................................Alternate Allele frequency filter.............................................................................
pos_retained = filtering.altAlleleFreqFilter(df,lklhd_df,minAltAlleleFreq)

# updating the read count matrix
df = df[pos_retained,:]
df_copy = df_copy[pos_retained,:]
alt_alleles = np.array(alt_alleles)
alt_alleles = alt_alleles[pos_retained]
lklhd_df = lklhd_df[pos_retained,:]
if isQualityFilteringOn:
   gq = gq[pos_retained,:]


# .................................Variant Removal filter.............................................................................
pos_retained_after_third = filtering.variantRemovalFilter(df,varThreshold)

# updating the read count matrix
df = df[pos_retained_after_third,:]
df_copy = df_copy[pos_retained_after_third,:]
alt_alleles = np.array(alt_alleles)
alt_alleles = alt_alleles[pos_retained_after_third]
lklhd_df = lklhd_df[pos_retained_after_third,:]
if isQualityFilteringOn:
   gq = gq[pos_retained_after_third,:]




# .................................Variant Removal mutated filter.............................................................................
final_pos_retained = []
final_pos_retained = filtering.variantRemovalMutated(df,lklhd_df,minMutFraction,alt_alleles,final_pos_retained)

end_time = time.time()

# updating the read count matrix

df = df[final_pos_retained,:]
df_copy = df_copy[final_pos_retained,:]
alt_alleles = alt_alleles[final_pos_retained]
if isQualityFilteringOn:
   gq = gq[final_pos_retained,:]

print("-------Filters applied !-------")
print("Time taken for filtering: ",end_time-filters_start_time)
print("\n")



# .................................Calculating likelihood at the remaninig sites.............................................................................
pd.DataFrame(df).to_csv(outputPrefixName+'modified_df_during_filters.tsv',sep='\t',header=False,index=False)
df[:,3] = alt_alleles
lklhd_df = lklhd_df[final_pos_retained,:]
pd.DataFrame(df).to_csv(outputPrefixName+'final_df.tsv',sep='\t',header=False,index=False)
print('Computing lklhd after filters....')
for i in range(df.shape[0]):
    counts = getAltAlleleCounts(df[i][4:])
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)
    lklhd_df[i] = np.where(coverage==0,np.nan,lklhd_df[i])

print("-------Final dimensions after all filters are applied-------")
print("Read count matrix shape : ",df.shape)
if isQualityFilteringOn:
   print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# # # # # # # # # # # # # # # # # # # # # # GRAPH LAPLACIAN # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
lklhd_df = lklhd_df.T
data_df = pd.DataFrame(lklhd_df)
data_df =KNNImputer(n_neighbors=8).fit_transform(X=data_df)
pd.DataFrame(data_df).to_csv(outputPrefixName+'final_lklhds.tsv',sep='\t',header=False,index=False)




lklhd_computed = pd.DataFrame(data_df.T)
n1 = lklhd_computed.shape[1]
l1 = lklhd_computed.shape[0]
zero_Ls = np.zeros((l1,n1))
one_Ls = np.zeros((l1,n1))
(uL, one_Ls, zero_Ls) = maximum_likelihood_tree.ML_initialization(lc=lklhd_computed, one_Ls=one_Ls, zero_Ls=zero_Ls)
lklhd_s = []
best_mat_s = []
label_s = []
inferredClusterLabelFileName = outputPrefixName+"inferred_cluster_labels"+".txt"
best_tree_s = []
best_cluster_s = []
ysc_s = []
sil_score_s=[]
best_geno_name_s = []

# # # # # # # # # # # # # # # # # # # # # # FINDING BEST LIKELIHOOD # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

if clustering == 'spectral':
    graph_laplacian = generate_graph_laplacian(df=data_df, nn=20)
    start = time.time()
    eigenvals, eigenvcts = compute_spectrum_graph_laplacian(graph_laplacian)
    print("Time taken for computing laplacian: ",time.time()-start," secs")


    eigenvcts_norms = np.apply_along_axis(
    lambda v: np.linalg.norm(v, ord=2), 
    axis=0, 
    arr=eigenvcts
    )

    eigenvals_sorted_indices = np.argsort(eigenvals)
    eigenvals_sorted = eigenvals[eigenvals_sorted_indices]
    index_lim = len(eigenvals_sorted[eigenvals_sorted<=eigenValueThreshold])

    # Eigen gap heuristic
    y=eigenvals_sorted[: index_lim]
    new_strategy = np.argmax(np.diff((y)))+1
    n_ev = new_strategy
    n_clusters = n_ev
    print('Number of candidate clusters after applying eigen gap heuristic: ',n_ev)

    clustering = SpectralClustering(n_clusters=n_clusters, eigen_solver = "arpack",
    affinity = 'rbf',assign_labels="discretize",random_state=42, n_init = 50).fit(data_df)
    labels = clustering.labels_
    final_labels = labels

    n_cluster = len(set(final_labels))
    print("\nNumber of clusters :",n_cluster)
    start_cluster = n_cluster - lt
    end_cluster = n_cluster + gt
    start_cluster = max(2,start_cluster)

    label_s.append(final_labels)
    for n_cluster in range(start_cluster, end_cluster+1):
        (best_L,best_mat,best_tree,best_geno_name) = maximum_likelihood_tree.findBestLikelihood(iterations,n_cluster,data_df,one_Ls,zero_Ls,final_labels,lklhd_computed)
        lklhd_s.append(best_L)
        best_mat_s.append(best_mat)
        best_tree_s.append(best_tree)
        best_geno_name_s.append(best_geno_name)

if clustering=='leiden':
    for resol in [0.1,0.2,0.3,0.4,0.5,0.75,1]:
        print("At resolution :",resol)
        adata = sc.AnnData(data_df)
        sc.pp.neighbors(adata,use_rep = 'X')
        sc.tl.leiden(adata,resolution=resol)
        cluster_labels = np.array(adata.obs['leiden'],dtype='int')
        final_labels=cluster_labels
        sil_score = silhouette_score(data_df, cluster_labels)
        label_s.append(cluster_labels)
        n_cluster = len(set(cluster_labels))
        print("Number of Clusters : ",n_cluster)
        if n_cluster in best_cluster_s:
            print("This cluster is already done")
            continue
        if n_cluster > 13:
            print("Number of clusters greater! skipping!")
            continue

        (best_L,best_mat,best_tree,best_geno_name) = maximum_likelihood_tree.findBestLikelihood(iterations,n_cluster,data_df,one_Ls,zero_Ls,final_labels,lklhd_computed)
        
        print("Resolution : ",resol,"\tLikelihood :",best_L,"\tSilhouette Score :",sil_score)
        best_cluster_s.append(n_cluster)
        lklhd_s.append(best_L)
        best_mat_s.append(best_mat)
        sil_score_s.append(sil_score)
        best_tree_s.append(best_tree)
        best_geno_name_s.append(best_geno_name)



print('Algorithm completed in: ',time.time()-start_time,' secs')
inferred_index = np.argmax(np.array(lklhd_s))
best_mat = best_mat_s[inferred_index]
best_labels = label_s[inferred_index]
final_inferred_clusters = len(set(best_labels))
best_tree = best_tree_s[inferred_index]
best_tree.write(format=9, outfile=outputPrefixName+"inferred_tree.nw")
best_geno = best_geno_name_s[inferred_index]


pd.DataFrame(best_geno).to_csv("Genotype configuration.tsv",sep='\t',index=False,header=False)
np.savetxt(inferredClusterLabelFileName,best_labels,fmt="%s", delimiter=",")
pd.DataFrame(best_mat).to_csv(inferredGenotypeFileName,index=False,header=False)
best_mat_trimmed = pd.DataFrame(best_mat).drop_duplicates()
inferred_genotypes = best_mat.T


final_indices = write_vcf.getRightOrderForVCF(df)
df = df[final_indices,:]
if isQualityFilteringOn:
   gq = gq[final_indices,:]
inferred_genotypes = inferred_genotypes[final_indices,:] 
vcfResultFileName = outputPrefixName+'outputInference.vcf'
vcfFile = open(vcfResultFileName,'w')
numCells = df.shape[1]-4
positions = df.shape[0]


print("------- SUMMARY -------")
print('Number of cells : ',df.shape[1]-4)
print('Number of genomic positions retained : ',df.shape[0])
print("Read count matrix shape : ",df.shape)
if isQualityFilteringOn:
   print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Inferred genotypes matrix : ",inferred_genotypes.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("-------------------------------------------------------")
print("-------------------------------------------------------")

print('Writing inferred genotypes as VCF file...')

if isQualityFilteringOn:
   write_vcf.writeVCFFileGQ(df,gq,vcfFile,inferred_genotypes)
else:
   write_vcf.writeVCFFile(df,vcfFile,inferred_genotypes)

print('Final inferred clusters in the tree: ',final_inferred_clusters)
print("Final inferred distinct genotypes: ",best_mat_trimmed.shape[0])

def generate_seaborn_colors(n): 
    palette = sns.color_palette("husl", n) 
    rgb_colors = [(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in palette] 
    return rgb_colors 
 
# Example usage 


data_df =pd.read_csv(outputPrefixName+'final_lklhds.tsv',sep='\t',header=None)
cluster = np.loadtxt(inferredClusterLabelFileName)
cluster= cluster.astype(int)
adata = sc.AnnData(X=data_df)
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X')
sc.tl.umap(adata,   min_dist = min_dist_defined)

n= len(set(cluster))
distinct_colors = generate_seaborn_colors(n) 



random_colors=dict()
for i in range(len(set(cluster))):
   random_colors[i] = distinct_colors[i]
colors=dict()
for i in range(len(cluster)):
   colors[i] = random_colors[cluster[i]]
labels = pd.DataFrame(colors.items(), columns=['index', 'label'])
adata.obs['batch'] = cluster.astype('str')
adata.uns['batch_colors'] = random_colors.values()
sc.pl.umap(adata,color='batch',show=False)
plt.savefig(outputPrefixName+"_umap.pdf")




with open('Final clusters inferred.txt','w') as f:
   f.write('Final inferred clusters in the tree\tFinal inferred distinct genotypes')
   f.write('\n')
   f.write(str(final_inferred_clusters) + '\t' + str(best_mat_trimmed.shape[0]))

print("Making GV tree...")
os.system("python3 GV_tree.py -out "+outputPrefixName)
