## Running PHALCON on simulated data

To generate simulated data, refer to ([Generate synthetic datasets](#generate-synthetic-datasets)).

### Parameters
#### Required parameters

`-i` : Input read count file name (*.tsv file*)

* PHALCON accepts data in the form of a readcounts file. Each row corresponds to a genomic site. The first four columns respectively correspond to *chromosome*, *genomic site*, *reference nucleotide*, and *alternate nucleotide*. The rest of the columns are the cells. Each entry is a 4-tuple format corresponding to the number of reads containing A,C,G,T (respectively) that aligned with that particular locus. *The file must be tab-separated.* (<font color="red">Dimension : *sites* x *(cells+4)*</font>)


##### Input format example
```
chr1	1	G	T	0,0,15,0	0,0,12,0	0,0,37,0	0,0,32,0	0,0,13,0	0,0,29,0	0,0,19,0	0,0,22,0	0,0,16,1	0,0,31,0
chr1	2	G	C	0,0,14,0	0,0,13,0	0,0,41,0	0,0,33,0	0,0,16,0	0,0,26,0	0,0,17,0	0,0,20,0	0,0,21,0	0,0,28,0
chr1	3	C	G	0,14,0,0	0,13,0,0	0,27,0,0	0,24,0,0	0,17,0,0	0,28,0,0	0,20,0,0	0,20,0,0	0,22,0,0	0,27,0,0
chr1	4	G	A	0,0,16,0	0,0,13,0	0,0,31,0	0,0,33,0	0,0,17,0	0,0,27,0	0,0,21,0	0,0,25,0	0,0,20,0	0,0,34,0
chr1	5	G	A	0,0,12,0	0,0,14,0	0,0,25,0	0,0,31,0	0,0,16,0	0,0,25,0	0,0,21,0	0,0,23,0	0,0,15,0	0,0,30,0
chr1	6	G	A	0,0,15,0	0,0,14,0	0,0,25,0	0,0,30,0	0,0,15,0	0,0,25,0	0,0,23,0	0,0,21,0	0,0,21,0	0,0,26,0
chr1	7	G	C	0,0,14,0	0,0,13,0	0,0,23,0	0,0,31,0	0,0,17,0	0,0,25,0	0,0,15,0	0,0,26,0	0,0,22,0	0,0,31,0
chr1	8	T	G	0,0,0,16	0,0,0,12	0,0,0,26	0,0,0,33	0,0,0,11	0,0,0,26	0,0,0,21	0,0,0,21	0,0,0,20	0,0,0,32
chr1	9	G	A	1,0,13,0	0,0,13,0	1,0,26,0	0,0,28,0	0,0,17,0	0,0,28,0	0,0,19,0	0,0,29,0	0,0,21,0	0,0,31,0
chr1	10	G	T	0,0,16,0	0,0,13,0	0,0,22,0	0,0,34,0	0,0,16,0	0,0,30,0	0,0,23,0	0,0,23,0	0,0,19,0	0,0,32,0
chr1	11	G	T	0,0,15,0	0,0,10,0	0,0,31,0	1,0,32,0	0,0,16,0	0,0,23,0	0,0,18,0	0,0,24,0	0,0,21,0	0,0,33,0
```

( **NOTE** : Panel seqeuncing outputs are generally in Loom format. Use ```loomToReadcount.py``` present in the ```supplementary``` folder to convert loom file to PHALCON-based format. To include indels, use ```loompyToReadcount_indels.py```)

#### Output parameters
`-o` : Output prefix (*string*)

* The output prefix for all the output files you will obtain after running PHALCON.

#### Optional parameters
`-gq` : Enable genotype quality filter (Default : 0) (*Optional*)

* If you have a genotype quality file and you want to use that for quality filtering, add ```-gq 1``` in the command while running PHALCON.

`-g` : Genotype quality file name (*.tsv file*) (*Optional*)

*  To use the genotype quality filtering, you need to pass the genotype quality file. Also, don't forget to enable the ```-gq``` parameter to not run into bugs.

The genotype quality file should look like:
```
0	75	14	21	33	86	18	7	71	93	20	78	23	50	84	53	55	94	95
1	84	84	65	89	78	8	67	12	6	13	84	95	87	23	89	47	29	49
2	59	19	45	24	9	84	71	61	26	87	67	7	11	14	7	11	11	25
3	32	94	17	90	47	28	9	37	26	80	6	81	18	24	49	48	22	93
4	25	38	58	25	88	82	25	41	15	49	81	78	41	80	77	53	57	80
5	30	88	83	84	76	93	39	57	34	87	54	33	63	88	58	9	19	40
6	57	27	95	96	77	51	51	49	43	39	30	18	6	94	14	62	53	73
7	25	10	78	55	23	36	48	54	96	44	41	6	49	44	7	18	43	32
8	83	65	33	98	7	68	94	44	77	72	69	48	30	81	50	30	47	25
9	13	91	32	20	52	13	28	31	39	21	60	85	9	56	55	41	6	76
10	62	80	64	44	72	71	43	60	70	40	84	95	25	43	92	90	87	67
11	80	85	10	54	61	20	69	35	98	39	47	81	43	44	69	59	6	71
12	17	91	95	96	32	41	78	67	16	56	5	5	95	35	71	57	59	60
```

The first column is the index, the rest of the columns mention the quality of the mapped reads. It is often present in the loom file. Extracting the genotype quality information for each site and cell is also embedded in the script ```loomToReadcount.py```. (<font color="red">Dimension:*sites* x *cells*</font>)



### Installation
Clone the repo and navigate to the main directory:
```bash
git clone https://github.com/Zafar-Lab/PHALCON
cd PHALCON
```

Create and activate the PHALCON conda environment:
```bash
conda env create -f environment.yml
conda activate phalcon
```
Verify that PHALCON and all dependencies have been installed successfully, run:
```python
python src/phalcon.py --help
```
This should display the list of available command-line arguments.

###Run PHALCON
The example below runs PHALCON using both a read count matrix (*sample_read_count_file.tsv*) and a genotype quality matrix (*sample_geno_qual_file*) while keeping all other parameters at their default values.

```python
python src/phalcon.py -i ../sample_read_count_file.tsv -g ../sample_geno_qual_file.tsv -gq 1
```

###Output

PHALCON mainly outputs the variant calls on each cell (.vcf format) and the reconstructed phylogeny (.gv format). Other auxiliary files, such as umap, cluster labels, etc, are also outputted.


### PHALCON executable
You may also create a system-wide executable for PHALCON. This allows PHALCON to be invoked directly from the command line using the `phalcon` command instead of `python src/phalcon.py`. To create the executable, navigate to the `src` directory and follow the steps below:

**Step 1**- Change directory to src folder and run the following command:
```bash
chmod +x phalcon.py
```

**Step 2**- Convert phalcon into an executable file using the following commands (second command is optional):
```bash
sudo mv phalcon.py /usr/local/bin/phalcon
sudo ln -s $(pwd)/phalcon.py /usr/local/bin/phalcon
```

After creating the symbolic link, PHALCON can be run directly as:
```bash
phalcon -i sample_read_count_file.tsv -g sample_geno_qual_file.tsv -gq 1
```

### Guidelines
#### Guidelines for filtering

The paramters related to filtering thresholds are described in detail below. **The defaults follow recommended best practices for such platforms, so they should not be changed unless necessary.**
 
`-r` : Minimum read depth threshold (*Default : 5*)

* Removes read count information from the cells with coverage depth less than this given threshold. For high coverage datasets, you can increase the value for more confident calls.

`-gq` and `-g` : Genotype quality filter (*Optional*, *Default Genotype quality : 30*)

* Removes read count information from cells with genotype quality less than 30. (Refer to [Optional parameters](#optional-parameters))


`-a` : Alternate frequency threshold (*Default : 0.2*)

* Filters out read-count information based on low variant allele frequency (VAF), since a VAF value much lower than 0.5 (for a heterozygous variant, expected VAF ~ 0.5) can indicate a potential false positive.

`-v` : Threshold for proportion of cells with insufficient read count information (*Default : 0.5*)

* Removes genomic sites where read count information is available for fewer than the given threshold proportion of cells.

`-m` : Threshold for proportion of sites harboring a mutation (*Default : 0.004*)

* Removes a genomic site if the fraction of cells harboring a mutation is very low. For real datasets, a higher cutoff (e.g., 0.01) may be used for more stringent filtering.

**A detailed analysis for parameter sensitivity related to some of these filters is present in the Supplementary document of the manuscript.**

#### Guidelines for clustering

`-c` : Clustering technique (*Default : spectral*)

* We specifically utilize graph-based clustering as it does not pre-assume any cluster structure, making it better suited for capturing heterogeneity. 

* You can choose either Spectral clustering or Leiden clustering to cluster the cells. Use: ```-c leiden``` for enabling leiden clustering.

### Generate synthetic datasets

To generate simulated data, use the script present in the folder [generate_synthetic_data](https://drive.google.com/drive/folders/1JWyuiIxNzRqUMV4MdHjfIVD-Hcdbi8ka?usp=drive_link).

You can provide the following arguments:

`--numCells` : Number of cells _(Default: 2000)_

* Specifies the total number of single cells to simulate.

`--numClones` : Number of clones _(Default: 10)_

* Specifies the number of distinct clonal populations.

`--numPos` : Number of genomic positions _(Default: 40000)_

* Specifies the total number of genomic positions to simulate.

`--covMean` : Mean sequencing coverage _(Default: 25)_

* Sets the average sequencing depth.

`--covVar` : Coverage variance _(Default: 50)_

* Controls the variability in sequencing coverage across sites.

`--dropoutRate` : Allelic dropout rate _(Default: 0.2)_

* Specifies the probability of allelic dropout during amplification.

`--mdaErrorRate` : Amplification error rate _(Default: 0.0025)_

* Specifies the error rate introduced during multiple displacement amplification.

`--sequencingErrorRate` : Sequencing error rate _(Default: 0.001)_

* Specifies the probability of sequencing errors.

`--seed` : Random seed _(Default: 129)_

`--prefixName` : Output file prefix _(Default: sc_2000_)_

* Prefix used for naming all generated output files.

**Note** : The dataset generator provides several additional parameters (e.g., average region length, amplification coefficients, and other simulation controls). For a complete list of available options, run:
```
python simulate_sc_datasets.py -h
```
