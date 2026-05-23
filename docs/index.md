# PHALCON

PHALCON is a scalable single-cell variant caller designed for high-throughput sequencing data. It is robust to common single-cell sequencing (SCS) errors and enables accurate mutation detection across large numbers of cells within practical runtimes.

![PHALCON overview](images/Overview_PHALCON.png)


## AML-67-001 tutorial

### Prepare the parameters
* `-D` : Tag for the dataset (*OPTIONAL*)(*Acts as an identifier*)
* `-i` : Input read count file name
* `-g` : Quality file name (*OPTIONAL*)
* `-t` : Enable/Disable Genotype Quality filtering (*default True*)
* `-o` : Output prefix

```python
data = 'AML_67_001'  # dataset
read_count_file_name = "/content/drive/MyDrive/PHALCON_git/aml_sample_run/Count_file_and_gq_file/ReadCounts_"+data+".tsv"  # readcount file
quality_file_name = "/content/drive/MyDrive/PHALCON_git/aml_sample_run/Count_file_and_gq_file/GQ_"+data+".tsv"   # genotype quality file
outputPrefix = '/content/drive/MyDrive/PHALCON_git/aml_sample_run/output/'+data+'_'    # output prefix
```

#### Additional parameters you can toggle with:

* `-r` : Minimum read deph for Read depth filter *(default 10)*
* `-q` : Minimum read quality threshold for Quality depth filter *(default 30)*
* `-a` : Minimum VAF for Alternate Allele Freq filter *(default 0.2)*
* `-v` : Minimum fraction of cells with available read count information required to retain a variant *(default 0.5)*
* `-m` : Minimum fraction of mutated cells required for a site to be retain a variant *(default 0.01)*
* `-s` : Seed

### Running PHALCON

 * The output obtained from MissioBio Tapestri pipeline is in the form of loom files. The loom file for AML-67-001 is present in the [AML-67-001 tutorial drive link](https://drive.google.com/drive/folders/1HAozFYcfdEZU5qF4-gxHTCm8h9Nl4hGI?usp=sharing) under the folder > "loom file"
* To convert the loom file in PHALCON-style input use the script ```loompy_to_readcount_indels.py``` present in the same folder

```python
!time python "/content/drive/MyDrive/PHALCON_git/aml_sample_run/python scripts for running phalcon/Algorithm_leiden_silhouette_on_indels.py"\
 -D $data -i $read_count_file_name -g $quality_file_name -o $outputPrefix
```

### Post-processing using dbSNP

The post processing includes removing variants with no mutation across all cells and removing clonal variants present in the dbSNP database 📑 

*  The dbSNP.vcf is present at the link [dbSNP_database](https://drive.google.com/file/d/1yy30skLXLOd4jDniXLWgaEqqozATlCN0/view?usp=drive_link)

Since the dbSNP data file is large, running it in COLAB kills the kernel.

Follow the steps below to generate the final vcf file of post-processed variants:

*  Download the folder [post-processing](https://drive.google.com/drive/folders/1pEHDfl8mCl3LuR02_PLioJCtQRRQ97Xx?usp=drive_link) in your local system
*   Open termnal inside the folder and type ```jupyter-notebook```
* Open the notebook [Post-processing_AML-67-001.ipynb](https://drive.google.com/file/d/1Yt3M_x0wBndzOCs9wl16o0pzKw_LAxfu/view?usp=drive_link)
* Run ```dbsnp_for_aml_finite_sites.py```
* You will get a final list of variants in form of a vcf file in a file name ending with ```_post_processed_final.vcf```

**NOTE**: As input to the ```dbsnp_for_aml_finite_sites.py``` python file, you only need the vcf file you obtained after running PHALCON. 
If required, change the location of the vcf file accordingly. 

* Use ```-i``` to provide location of the input vcf file 
* Use ```-o``` to mention the name of the post processed vcf file

### Annotate variants

For variant annotation using ANNOVAR, the necessary files and scripts are present in the [annotating_variants](https://drive.google.com/drive/folders/1pfyQ7VWEw9cZKXe9RLwVetxKC4SY2vss?usp=drive_link) link. 


Since running ANNOVAR in COLAB kills the kernel, use the following steps to annotate mutations:

*  Download the folder [annotating_variants](https://drive.google.com/drive/folders/1pfyQ7VWEw9cZKXe9RLwVetxKC4SY2vss?usp=drive_link) in your local system
*   Open termnal inside the folder and type ```jupyter-notebook```
* Open the notebook [Annotating_the_final_variants.ipynb](https://drive.google.com/file/d/1WBz1aBEFTjcZbZ4y6XOBEubgQwR3gvGk/view?usp=drive_link)
* Run ```run_annovar_for_annotating_final_variants.py```
* You will get an annotated files of the final variants

**NOTE**: As input to the ```run_annovar_for_annotating_final_variants.py``` python file, you only need the post-processed vcf file. 
If required, change the location of the vcf file accordingly. 

* Use ```-i``` to provide location of the post-processed input vcf file 
* Use ```-o``` to mention the output tag of the annotated files









