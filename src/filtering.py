# IMPORTS ###########################################################################################################
import argparse
import numpy as np
from tqdm import tqdm
#####################################################################################################################


argParser = argparse.ArgumentParser(prog='PROG')

argParser.add_argument('-s', '--seed', type=int, default = 12099)



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


def genoQualityFilter(df,gq,qualityVal):
# Genotype quality filter 
# Puts the readcount as "0,0,0,0" at places where the genotype quality is lower than a threshold
    muts = df.shape[0]
    cells = df.shape[1]-4
    for i in tqdm(range(muts),desc="Applying GENO QUALITY filter..."):
        for j in range(4,4+cells):
            if gq[i][j-4] < qualityVal:
                if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                   df[i][j] = "0,0,0,0"
                elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                   df[i][j] = "0,0"
                else:
                   print("Something wrong in gq filter at pos",i," cell ",j)



def readDepthFilter(df,minReadDepth):
# Read depth quality filter
# Puts the readcount as "0,0,0,0" at places where the read depth quality is lower than a threshold
    muts = df.shape[0]
    cells = df.shape[1]-4
    for i in tqdm(range(muts),desc="Applying READ DEPTH QUALITY filter..."):
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        read_depth = np.sum(counts,axis=1)
        for j in range(4,4+cells):
            if read_depth[j-4] < minReadDepth:
                if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                   df[i][j] = "0,0,0,0"
                elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                   df[i][j] = "0,0"
                else:
                   print("Something wrong in read depth filter at pos",i," cell ",j)



def altAlleleFreqFilter(df,lklhd_df,minAltAlleleFreq):
# Alternate Allele frequency filter (To remove FP due to WGA)
# If there is no alternate allele, put "0,0,0,0" at all cells across the readcount dataframe
# Else, if likelihood > 0.5 and alternate/total depth is less than a certain threshold 
#       then put "0,0,0,0" at that place in the readcount dataframe
    muts = df.shape[0]
    cells = df.shape[1]-4
    pos_retained_list = []

    for i in tqdm(range(muts),desc="Applying ALTERNATE ALLELE FREQUENCY filter..."):
        if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
            allele_counts = np.zeros((cells,4),dtype=int)
            ref = df[i][2]
            ref_ind = charToIndex(ref)
            counts = getAltAlleleCounts(df[i][4:])
            counts = np.array(counts.tolist())
            allele_counts = np.sum(counts,axis=0)  # Total allele count of each base at a certain site. Shape : (4,)
            ref_count = allele_counts[ref_ind]
            allele_counts[ref_ind] = -1
            alt_allele = 'X'
            if np.max(allele_counts) > 0:
                alt_ind = np.argmax(allele_counts)
                alt_allele = indexToChar(alt_ind)
                pos_retained_list.append(i)
        elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
            allele_counts = np.zeros((cells,2),dtype=int)
            counts = getAltAlleleCounts(df[i][4:])
            counts = np.array(counts.tolist())

            allele_counts = np.sum(counts,axis=0)  # Total allele count of each base at a certain site. Shape : (2,)
          
            alt_allele = 'X'
            allele_counts[0] = -1  # in case of indels, reference index is always zero and alternate index is always 1 because its counts go like ref:alt
            if np.max(allele_counts) > 0:
                alt_ind = np.argmax(allele_counts)
                alt_allele = df[i][3]
                pos_retained_list.append(i)
        else:
            print("Something wrong in alt allele freq filter at pos ",i)

        if alt_allele=='X':      
            if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
               df[i][4:] = np.full(cells, '0,0,0,0')
               continue
            elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
               df[i][4:] = np.full(cells, '0,0')
               continue
            else:
               print("Something wrong with alternate allele at alt allele freq filter at pos ",i)
        else:
            for cell in range(cells):
                counts = np.array(list((map(int,df[i][4+cell].strip().split(",")))))
                if df[i][4+cell] == '0,0,0,0' or df[i][4+cell] == '0,0':
                    continue
                read_depth = sum(counts)
                if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                   alt_depth = counts[alt_ind]
                elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                   alt_depth = counts[1]  # because indels carry readcount info as ref:alt, so reference depth is present at index 0 and alternate depth is present at index 1   

                if lklhd_df[i][cell]>0.1 and np.float16(alt_depth/read_depth) < minAltAlleleFreq:
                    if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                       df[i][4+cell] = '0,0,0,0'
                    elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                       df[i][4+cell] = "0,0"
                    else:
                       print("Something wrong with read count information at later part of alt allele freq filter at pos ",i," and cell ",cell)
    return pos_retained_list
                    
                

def variantRemovalFilter(df,threshold):
# Variant removal filter (based on read count information across all cells)
# If read count is  not "0,0,0,0" across more than a certain threshold among all cells then retain that site
    muts = df.shape[0]
    cells = df.shape[1]-4
    pos_retained = []

    for i in tqdm(range(muts),desc="Applying VARIANT REMOVAL filter..."):
        cell_count = 0
        for cell in range(cells):
            if df[i][4+cell] not in ['0,0,0,0','0,0']:
                cell_count+=1
        if cell_count >= (threshold*cells):
            pos_retained.append(i)
    print('num of pos retained after third filter: ',len(pos_retained))
    return pos_retained



def variantRemovalMutated(df,lklhd_df,minMutFraction,alt_alleles,pos_retained):
# Variant removal filter (based on fraction of cells mutated)
# If the likelihood is greater than 0.5 and read count is NOT "0,0,0,0" then increase the cell count by 1
#    If cell count is greater than a certain threshold then retain that site
    muts = df.shape[0]
    cells = df.shape[1]-4


    for i in tqdm(range(df.shape[0]),desc="Applying VARIANT REMOVAL MUTATED filter..."):
        cell_count = 0
        for cell in range(cells):
            if lklhd_df[i][cell]>0.5 and df[i][4+cell] not in ['0,0,0,0','0,0']:
                cell_count+=1
        if cell_count >= (minMutFraction*cells):
            pos_retained.append(i)
    print("Length of positions retained after FINAL FILTER: ",len(pos_retained))
    for pos in pos_retained:
        print(df[pos][0],"\t",df[pos][1],"\t",df[pos][2],"\t",' alt inferred: ',alt_alleles[pos],' true_alt: ',df[pos][3])
    return pos_retained


getAltAlleleCounts =np.frompyfunc(getAltAlleleCounts, 1, 1)  
