# IMPORTS ###########################################################################################################
import sys
import numpy as np
import re
#####################################################################################################################




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



def getSortedChrPos( keys ):
# Sorting positions for sending it to the VCF file
    key_func = lambda txt: [charToNum(c) for c in re.split('([0-9]+)', txt)]
    return sorted(keys, key = key_func)



def getRightOrderForVCF(df):
# Returns the final list of indices to be sent to the vcf file
    final_indices = []
    unique_chr_keys = np.unique(df[:,0])
    chr_keys = np.array(getSortedChrPos(unique_chr_keys))
    dict={}
    for key in chr_keys:
        dict[key] = []
    for row in range(df.shape[0]):
        dict[df[row][0]].append(row)
    for key in dict.keys():
        for ind in dict[key]:
            final_indices.append(ind)
    return final_indices



def writeVCFHeader(vcfFile,numCells):
# Writes the header of the VCF file
    vcfFile.write("##fileformat=VCFv4.1\n")
    vcfFile.write("##source=OurAlgo" + "OurAlgo v" + '0' + "." + '1' + "." + '0' + "\n")
    vcfFile.write("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
    vcfFile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n")
    vcfFile.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for alt alleles\">\n")
    vcfFile.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
    vcfFile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
    vcfFile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcfFile.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
    vcfFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    
    for i in range(1,numCells+1):
        vcfFile.write("\tcell"+str(i))
    vcfFile.write('\n')



def writeEntryGQ(df,gq,vcfFile,inferred_genotypes):  
# df = Read count dataframe (sites x cells)
# Writes entry of the final positions in the VCF file
    pos = df.shape[0]
    numCells = df.shape[1]-4
    alt_depth = 0
    count_entry = []
    for i in range(pos):
        vcfFile.write(str(df[i][0])+'\t')  # Chr
        vcfFile.write(str(df[i][1])+'\t')  # Position
        vcfFile.write('*\t')
        vcfFile.write(str(df[i][2])+'\t')  # Reference allele
        vcfFile.write(str(df[i][3])+'\t')  # Alternate allele
        vcfFile.write('*\t')
        vcfFile.write('PASS\t')
        vcfFile.write('DP=')
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        allele_counts = np.sum(counts,axis=0)
        depthAllCells = np.sum(allele_counts)
        vcfFile.write(str(depthAllCells)+'\t')
        vcfFile.write('GT:AD:DP:GQ:PL\t')
        for cell in range(numCells):
            if int(inferred_genotypes[i][cell])==1:
                vcfFile.write('0/1:')
                qual = str(gq[i][cell])
            else:
                vcfFile.write('0/0:')
                qual = str(gq[i][cell])
            count_entry = list(map(int,df[i][4+cell].strip().split(",")))
            alt_depth = int(count_entry[charToIndex(df[i][3])])
            vcfFile.write(str(alt_depth)+':')
            total_depth = str(np.sum(np.array(count_entry)))
            vcfFile.write(total_depth+':'+qual+'\t')
        vcfFile.write('\n')



def writeEntry(df,vcfFile,inferred_genotypes):  
# df = Read count dataframe (sites x cells)
# Writes entry of the final positions in the VCF file
    pos = df.shape[0]
    numCells = df.shape[1]-4
    alt_depth = 0
    count_entry = []
    for i in range(pos):
        vcfFile.write(str(df[i][0])+'\t')  # Chr
        vcfFile.write(str(df[i][1])+'\t')  # Position
        vcfFile.write('*\t')
        vcfFile.write(str(df[i][2])+'\t')  # Reference allele
        vcfFile.write(str(df[i][3])+'\t')  # Alternate allele
        vcfFile.write('*\t')
        vcfFile.write('PASS\t')
        vcfFile.write('DP=')
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        allele_counts = np.sum(counts,axis=0)
        depthAllCells = np.sum(allele_counts)
        vcfFile.write(str(depthAllCells)+'\t')
        vcfFile.write('GT:AD:DP:GQ:PL\t')
        for cell in range(numCells):
            if int(inferred_genotypes[i][cell])==1:
                vcfFile.write('0/1:')
                qual = str(80)
            else:
                vcfFile.write('0/0:')
                qual = str(80)
            count_entry = list(map(int,df[i][4+cell].strip().split(",")))
            alt_depth = int(count_entry[charToIndex(df[i][3])])
            vcfFile.write(str(alt_depth)+':')
            total_depth = str(np.sum(np.array(count_entry)))
            vcfFile.write(total_depth+':'+qual+'\t')
        vcfFile.write('\n')



def writeVCFFileGQ(df,gq,vcfFile,inferred_genotypes):
# Writes the whole VCF file (header + final positions)
    numCells = df.shape[1]-4
    writeVCFHeader(vcfFile, numCells)
    writeEntryGQ(df,gq,vcfFile,inferred_genotypes)
   


def writeVCFFile(df,vcfFile,inferred_genotypes):
# Writes the whole VCF file (header + final positions)
    numCells = df.shape[1]-4
    writeVCFHeader(vcfFile,numCells)
    writeEntry(df,vcfFile,inferred_genotypes)




def getAltAlleleCounts(arr):
# Converts read count strings (e.g. "0,1,2,1") to numpy array (e.g. [0,1,2,1])
    return np.array(list(map(int,arr.strip().split(","))))


getAltAlleleCounts =np.frompyfunc(getAltAlleleCounts, 1, 1)  
