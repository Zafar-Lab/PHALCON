
## TN4 tutorial

### Prepare the parameters
* `-D` : TN4
* `-i` : readCounts_MBTN4_nodups.tsv (Input read count file name , present in the folder [Count_file_and_gq_file](https://drive.google.com/drive/folders/14__kHyT4DfLZl_JTa-NJVHbDUGbbJ8Ze?usp=drive_link))
* `-g` : gq_file_TN4.tsv (Quality file name, present in the folder [Count_file_and_gq_file](https://drive.google.com/drive/folders/14__kHyT4DfLZl_JTa-NJVHbDUGbbJ8Ze?usp=drive_link))
* `-o` : TN4_ (Output prefix)


### Running PHALCON

* The output obtained from MissioBio Tapestri pipeline is in the form of loom files. The loom file for TN4 is present in the folder [loom file](https://drive.google.com/drive/folders/1vW4XBpDR5H7cTjCcf_3dkokAOKc9Otcr?usp=drive_link).
* To convert the loom file in PHALCON-style input use the script ```loomToReadcount.py``` present in the same folder
* Scripts for running PHALCON are present at [python scripts for running phalcon](https://drive.google.com/drive/folders/1Evck59kj2s5RYVdNzXWCGzAUIREU1_5A?usp=drive_link)

<font color="red">Input</font> : Pass the readcount file, genotype quality file and the output prefix name.
Run the following command:
```python
python Algorithm_hierarchical.py -D TN4 -i readCounts_MBTN4_nodups.tsv -g gq_file_TN4.tsv -o TN4_
```
<font color="red">Output</font> : 
You will get the variants inferred by PHALCON in form of a VCF file.

Subset vcf output:
```
##fileformat=VCFv4.1
##source=OurAlgoOurAlgo v0.1.0
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  cell1   cell2   cell3   cell4   cell5   cell6   cell7   cell8   cell9   cell10  cell11>
chr1    28316313        *       C       A       *       PASS    DP=853946       GT:AD:DP:GQ:PL  0/0:0:20:60     0/1:10:10:30    0/1:39:47:57    0/1:387:388:99>
chr1    40430992        *       G       C       *       PASS    DP=559770       GT:AD:DP:GQ:PL  0/0:2:40:70     0/1:82:150:99   0/1:30:53:99    0/1:172:437:99>
chr1    45811136        *       C       T       *       PASS    DP=1282176      GT:AD:DP:GQ:PL  0/1:25:50:99    0/1:32:56:99    0/1:148:262:99  0/1:502:1054:9>
chr1    62740446        *       T       C       *       PASS    DP=1900868      GT:AD:DP:GQ:PL  0/1:32:81:99    0/1:75:189:99   0/1:80:170:99   0/1:400:1084:9>
chr1    62740449        *       T       C       *       PASS    DP=1909598      GT:AD:DP:GQ:PL  0/1:32:81:99    0/1:75:189:99   0/1:80:170:99   0/1:405:1084:9>
chr1    62971394        *       G       C       *       PASS    DP=3745080      GT:AD:DP:GQ:PL  0/0:2:179:99    0/1:79:155:99   0/1:219:553:99  0/1:484:984:99>
chr1    89734417        *       C       G       *       PASS    DP=1321902      GT:AD:DP:GQ:PL  0/0:2:278:99    0/1:120:238:99  0/1:78:198:99   0/1:131:446:99>
chr1    94502814        *       G       C       *       PASS    DP=557482       GT:AD:DP:GQ:PL  0/0:0:24:72     0/1:8:33:84     0/1:14:90:47    0/1:38:117:99 >
chr1    145015876       *       A       G       *       PASS    DP=1541480      GT:AD:DP:GQ:PL  0/0:2:168:99    0/0:19:132:75   0/0:4:205:99    0/0:6:1224:99 >
chr1    145112420       *       C       T       *       PASS    DP=1947138      GT:AD:DP:GQ:PL  0/1:62:156:99   0/1:27:119:99   0/1:101:367:99  0/1:481:1447:9>  
```

### Post-processing and annotation

The intermediate steps to obtain the final annotated list of variants are explained in detail below:

#### Post-processing using dbSNP and matched normal removal

For TNBC, the post processing includes removing variants with no mutation across all cells, removing clonal variants present in the dbSNP database 📑 and removing variants present in the matched normal sample.

*  The dbSNP.vcf is present at the link [dbSNP_database](https://drive.google.com/file/d/1yy30skLXLOd4jDniXLWgaEqqozATlCN0/view?usp=drive_link)
(Change the location of the dbsnp vcf file accordingly in the python script.)

* Python script to perform post-processing is present in the same folder under the name `dbsnp_normal_removal_from_tnbc.py`

**NOTE**: As input to the `dbsnp_normal_removal_from_tnbc.py` python file, you only need the vcf file you obtained after running PHALCON. 
If required, change the location of the vcf file accordingly. 

* Use ```-i``` to provide location of the input vcf file 
* Use ```-o``` to mention the name of the post processed vcf file

#### Annotate variants

For variant annotation using ANNOVAR, the necessary files and scripts are present in the [annotating_variants](https://drive.google.com/drive/folders/12Kmuvd4ZcNbwtCSbGjrhKp4OudKsaZQ1?usp=drive_link) link. 

* Python script to obtain the annotations is present in the same folder under the name `run_annovar_for_annotating_final_variants.py`

**NOTE**: As input to the `run_annovar_for_annotating_final_variants.py` python file, you only need the post-processed vcf file. 
If required, change the location of the vcf file accordingly. 

* Use ```-i``` to provide location of the post-processed input vcf file 
* Use ```-o``` to mention the output tag of the annotated files

### Visualization

The genotyped variants are shown in the form of a heatmap, and the chronological order of mutations is shown in the form of a phylogenetic tree.

The following files are needed to produce the final tree and genotypes (All of these files can be obtained by running the above pipelines):

* ```-d``` : Tag for the dataset (*OPTIONAL*)(*Acts as an identifier*)
* ```--cluster_labels_file``` : Inferred cluster labels (*Obtained as one of the files after running PHALCON*)
* ```--lklhd_file``` : Likelihood file containing mutation likelihood for each cell and variant (*Obtained as one of the files after running PHALCON*)
* ```--final_df_file``` : Likelihood file containing mutation likelihood for each cell and variant (*Obtained as one of the files after running PHALCON*)
* ```--genotype_config_file``` : Genotype configuration file (*Obtained as one of the files after running PHALCON*)
* ```--newick_file``` : Newick format file (*Obtained as one of the files after running PHALCON*)
* ```--vcf_file``` : Final post-processed VCF file (Post-processed final VCF file, refer to [Post-processing using dbSNP](#post-processing-using-dbsnp))
* ```--variant_function_file``` : Gene-annotated file of the post-processed VCF file, refer to [Annotate variants](#annotate-variants))

**For a consolidated run including post processing, annotation and visualisation, follow the steps**:

* Download the folder [post_processing_plotting_and_visualisation](https://drive.google.com/drive/folders/1WEFGqmpW1zGKFy9tHxB_Olx5HT8XvVqg?usp=drive_link) in your local system.
* Open terminal inside the folder and type ```jupyter-notebook```.
* Open the notebook [Post_processing_anotation_and_visualisation](https://drive.google.com/file/d/17G8aQVumRexSpp_QBzuPTWnTYcELdC9K/view?usp=drive_link) and run the contents.
* The outputs will be saved in the [Output](https://drive.google.com/drive/folders/1dUHPSmZvjNXg105h9JqB6aSgEUeaJYYr?usp=drive_link) folder.

You will obtain the following outputs:

* <font color="red">**Genotype heatmap**</font>:

The genotype heatmap shows the genotype profiles of each cell. There is a bar at the top which shows the cluster association of each cell.
![Heatmap](images/TN4_heatmap.png)

* <font color="red">**UMAP**</font>:

The umap shows the clusters obtained by PHALCON on a 2-d plane.
<center>
<img src="../images/TN4_umap.svg" width="80%">
</center>

* <font color="red">**Phylogenetic tree**</font>:

The reconstructed mutation history tree shows the order in which the mutations have appeared.
![Phylogenetic_tree](images/TN4_clonal_tree.png)

