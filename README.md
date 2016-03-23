## Overview ##

**Monovar** is a single nucleotide variant (SNV) detection and genotyping algorithm for single-cell DNA sequencing data. It takes a list of bam files as input and outputs a vcf file containing the detected SNVs.

## Dependencies ##

* Python: NumPy v1.8.1 ([http://www.numpy.org/]()), SciPy v0.14.0 ([http://www.scipy.org/]()), Pysam v0.8.1 ([https://code.google.com/p/pysam/]())
* Samtools v0.1.19 (included in the external folder)

## Installation ##

Clone the Monovar repository: 

```
#!python
git clone git@bitbucket.org:hamimzafar/monovar.git
cd monovar

```
Install the Monovar python package:

```
#!python

sudo python setup.py install
```

Give execute permission to the file monovar.py

```
#!python

chmod +x src/monovar.py
```

Add the samtools folder and src folder to the PATH:

```
#!python

CURR_DIR=$(pwd)
export PATH=$PATH:$CURR_DIR/external/samtools
export PATH=$PATH:$CURR_DIR/src
```

## Usage ##
The program requires multiple bam files. The bam files should be sorted by coordinates. The raw sequence reads in .fastq format should be aligned to a reference genome with the help of an aligner program (e.g., BWA ([http://bio-bwa.sourceforge.net/]())). Aligner like BWA generates sam files containing aligned reads. The sam files can be converted to compressed bam files using ```samtools view``` command (see Samtools manual for details [http://www.htslib.org/doc/samtools.html]()). 

We have included three sample bam files in the folder examples. To run Monovar, a reference genome file is also needed. Assuming indexed reference genome file to be ref.fa and present in the examples directory, go to the examples directory and run Monovar on the provided bam files as follows:

```
#!python

samtools mpileup -BQ0 -d10000 -f ref.fa -q 40 -b filenames.txt | monovar.py -p 0.002 -a 0.2 -t 0.05 -m 2 -f ref.fa -b filenames.txt -o output.vcf
```
The arguments of Monovar are as follows:

```
#!python

-b: Text file containing the full path for each Bam file. One file per line.
-f: Reference genome file.
-o: Output file.
-t: Threshold to be used for variant calling (Recommended value: 0.05)
-p: Offset for prior probability for false-positive error (Recommended value: 0.002)
-a: Offset for prior probability for allelic drop out (Default value: 0.2)
-m: Number of threads to use in multiprocessing (Default value: 1)
-c: Flag indicating whether to use Consensus Filter (CF) or not (Possible values: 0, 1; Default Value: 1; if 1 then CF is used, otherwise not used)  
```
We recommend using cutoff 40 for mapping quality when using ```samtools mpileup```. To use the probabilistic realignment for the computation of Base Alignment Quality, drop the ```-B``` while running ```samtools mpileup```.