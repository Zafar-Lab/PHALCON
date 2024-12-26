# PHALCON
## Description
PHALCON is a variant calling algorithm specifically designed for high-throughput single-cell sequencing datasets.
## Usage
PHALCON takes as input a read count matrix $(sites \times cells)$. If you have a loom file instead, a script named ```loomToReadcount.py``` is present in the main directory, output of which can be fed as input to PHALCON.
For running PHALCON, download the ```src``` folder, unzip it on your system and follow the steps below:

**Step 1** Change directory to ```src``` folder and run the following command:
```
chmod +x phalcon.py
```
**Step 2** Convert phalcon into an executable file using the following commands (_second command is optional_):
```
sudo mv phalcon.py /usr/local/bin/phalcon
sudo ln -s $(pwd)/phalcon.py /usr/local/bin/phalcon
```
###Run phalcon by giving inputs on the command line
Below is an example where "sample_read_count_file.tsv" and "sample_gq_file.tsv" files are provided as input with all other variables being kept at the default values.
```
phalcon -i ../sample_read_count_file.tsv -g ../sample_gq_file.tsv
```
## Help
Run ```phalcon -help``` for the description of parameters, along with their default values.
