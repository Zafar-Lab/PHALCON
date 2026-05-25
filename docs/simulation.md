## Running PHALCON on simulated data

To generate simulated data, refer to refer to [Generate synthetic datasets](#generate_synthetic_datasets)).

### Know the parameters

Each parameter is explained in detail below:

* `-i` : Input read count file name (*.tsv file*)

PHALCON accepts data in the form of a readcounts file. The readcount file should look like:
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

Each row corresponds to a genomic site. The first four columns respectively correspond to *chromosome*, *genomic site*, *reference nucleotide*, and *alternate nucleotide*. The rest of the columns are the cells. Each entry is a 4-tuple format corresponding to the number of reads containg A,C,G,T (respectively) that aligned with that particular locus. (<font color="red">Dimension:*sites* x *cells*</font>)

( **NOTE** : Output obtained by panel seqeuncing datasets is generally in the form of a loom file. Use ```loomToReadcount.py``` present in the ```supplementary``` folder to convert loom file to PHALCON-based format. To include indels, use ```loompyToReadcount_indels.py```)

* ```-gq``` : Enable genotype quality filter (Default : 0) (*Optional*)

If you have a genotype quality file and you want to use that for quality filtering, add ```-gq 1``` in the command while running PHALCON.

* ```-g``` : Genotype quality file name (*.tsv file*) (*Optional*)

To use the genotype quality filtering, you need to pass the genotype quality file. Also, don't forget to enable the ```-gq``` parameter to not run into bugs.

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

* `-o` : Output prefix (*string*)
The output prefix for all the output files you will obtain after running PHALCON.
