# SVRE
Structural Variation detection using Relative Entropy

SVRE is a program that predicts structural variations (deletions, duplications, inversions, etc.) in a genome. One of the nice things about SVRE is that it has relatively low background noise levels, and thus predictions are relatively easy to identify visually from the visualizations that are created.

SVRE is designed to work with paired-end short read sequencing data. The algorithm behind how SVRE works is based on information theory and relative entropy.

# I'm in a hurry
Make sure you have sorted and indexed your bam files.
```
svre.pl -r1 <R1 bam file> -r2 <R2 bam file> -ori FR -output <output_prefix>
```
Then take a look at the `output_prefix_graph.png` image.

# Prerequisites
The following perl modules are required:
* Cwd
* File::Temp
* Getopt::Long
* List::Util
* Math::BigFloat
* Math::CDF
* Math::Round
* Math::Trig

The following additional software is required (and to be present on your default path):
* samtools
* Rscript

# Installation
There are only two pieces for this program - the main `svre.pl` perl script, and the `sv.pm` perl module. The script only needs to be executable (typically `chmod +x` on a \*nix system). `sv.pm` must be available - usually in linux this can be done by making sure it's in a directory that is found in your `PERL5LIB` environment variable (or adding the path to the folder containing `sv.pm` to that variable).

# Installation Recipe
This should get you up and running on a base Ubuntu 18.04 (tested on AWS AMI `ami-0dad20bd1b9c8c004` (Canonical, Ubuntu, 18.04 LTS, amd64 bionic image build on 2019-02-12)).
```
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install make gcc bwa samtools r-base-core libmath-round-perl
sudo cpan -i Math::CDF
cd /usr/local/src
sudo git clone https://github.com/swainechen/svre
sudo mkdir /usr/local/lib/site_perl
sudo ln -s /usr/local/src/svre/sv.pm /usr/local/lib/site_perl
sudo chmod +x /usr/local/src/svre/svre.pl
sudo ln -s /usr/local/src/svre/svre.pl /usr/local/bin
```
Then the example below should run (you need your own reference genome and fastq files).

# Basic background
SVRE requires two bam files as input.
These bam files should be created with either `bwa` or `bowtie` (version 1 or 2), as some knowledge of the mapping flags is required. We recommend `bwa` 0.7.10 or above.

Because most mappers actually use mapping distances to choose between possibilities for mapping position (i.e. they make assumptions about the distribution of the mapping distances), SVRE requires that each read gets mapped singly. Therefore, you need to map the R1 reads by themselves, and then separately map the R2 reads by themselves. These files are then the input to SVRE.

SVRE is also designed as a method that has good sensitivity, but this relies on having a good reference genome. Using a different reference, even from the same species, can degrade the signal to noise ratio dramatically.

As read length has increased with advances in second-generation sequencing technology, sequencing accuracy still often is lower towards the end of a read than at the beginning. To maintain sensitivity, SVRE only uses reads that map with 100% identity. This can filter out many reads if your read length is quite long. For most bacterial genomes, a good balance therefore is to trim your reads to the first 75 bp (or to trim by quality to this length). Then perform the mapping.

# Example
Let's say we have all this in a single directory:
* `reference.fna` - a fasta file of the reference sequence
* `R1.fastq` - first read of paired end sequencing data - trimmed to 75bp
* `R2.fastq` - second read of paired end sequencing data - trimmed to 75bp

Then we can do the following:
```
bwa index reference.fna
bwa aln reference.fna R1.fastq > R1.sai
bwa aln reference.fna R2.fastq > R2.sai
bwa samse reference.fna reference.fna R1.sai R1.fastq | samtools view -bS - > R1.bam
bwa samse reference.fna reference.fna R2.sai R2.fastq | samtools view -bS - > R2.bam
samtools sort R1.bam R1-sort
samtools sort R2.bam R2-sort
samtools index R1-sort.bam
samtools index R2-sort.bam
svre.pl -r1 R1-sort.bam -r2 R2-sort.bam -ori FR -output svre_results
```

# Output
Three files are output by the svre.pl program (svre_results is used here based on the example above, set this with the -output parameter):
* `svre_results_graph.png` - This gives an easy to view image of the relative entropies across all windows. Spikes in the relative entropy graph correspond with predicted structural variations.
* `svre_results_re.txt` - This gives summaries of all the windows and relative entropies.
* `svre_results_detail.txt` - This gives details of all the windows that were deemed to be significant.
* `svre_results_SV.txt` - This gives summarized predictions of structural variations detected. This is still very rough at this point.
