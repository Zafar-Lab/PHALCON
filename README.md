# PHALCON
## Description
PHALCON is a fast, efficient variant caller that is robust to SCS-errors. Specifically designed for high-throughput data, PHALCON performs reliable variant calling on a large number of cells within a reasonable timeframe. 

## Usage
PHALCON takes as input a read count matrix $(sites \times cells)$ and (_optionally_) a genotype quality matrix. If you have a loom file instead, a script named ```loomToReadcount.py``` is present in the main directory, output files of which can be fed as input to PHALCON.
For running PHALCON, download the ```src``` folder, unzip it on your system, and follow the steps below:

**Step 1**- Change directory to ```src``` folder and run the following command:
```
chmod +x phalcon.py
```
**Step 2**- Convert phalcon into an executable file using the following commands (_second command is optional_):
```
sudo mv phalcon.py /usr/local/bin/phalcon
sudo ln -s $(pwd)/phalcon.py /usr/local/bin/phalcon
```
## Arguments
```-i```: Input read count file

```-o```: Output prefix

```-r```: Minimum read depth threshold (_Default_ : 5)

```-a```: Alternate frequency threshold (_Default_ : 0.2)

```-v```: Threshold for proportion of cells with insufficient read count information (_Default_ : 0.5)

```-m```: Threshold for proportion of sites harboring a mutation (_Default_ : 0.004)


## Optional Arguments
```-gq```: Enable genotype quality filter

```-q```: Genotype quality threshold (_Default_ : 0.2)


## Run phalcon

On the command line, give the input arguments (use ```help``` for the list of arguments) and run phalcon.

Below is an example where "sample_read_count_file.tsv" and "sample_geno_qual_file.tsv" files are provided as input with all other variables being kept at the default values.
```
phalcon -i ../sample_read_count_file.tsv -g ../sample_geno_qual_file.tsv
```
Use ```-gq 0``` to disable the genotype quality filter. For a sample run, you can find the input files here:

[Sample read count and genotype quality files](https://drive.google.com/drive/u/1/folders/1DuhxBdxZNmsljerC1NVS12M_1r0B4Sbw)

## Help
Run ```phalcon -help``` for the description of parameters, along with their default values.
