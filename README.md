<p align="center">
  <img src="docs/chivariariba.png" alt="ChiVariARIBA Logo" width="600"/>
</p>

<p align="center">
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/language-Python-blue.svg" alt="Python Badge"/></a>
  <a href="https://www.gnu.org/software/bash/"><img src="https://img.shields.io/badge/shell-Bash-lightgrey.svg" alt="Bash Badge"/></a>
  <a href="https://creativecommons.org/licenses/by/4.0/"><img src="https://img.shields.io/badge/license-CC--BY-blue.svg" alt="License: CC-BY"/></a>
</p>

# ChiVariARIBA

Chitin Metabolism Pathway Gene Identification By Assembly Using ARIBA

This resource was developed to allow the user to effectively identify chitin metabolism pathway components in Vibrios and other related species using paired-end sequencing reads as input.

For how to use ARIBA, please see the [ARIBA wiki page][ARIBA wiki]. *THIS README.md FILE HAS BEEN ADAPTED WITH MODIFICATIONS FROM THE ARIBA README.md FILE*

## Contents
* [Introduction](#introduction)
* [Preparing ChiVariARIBA reference sequences](#Preparing-ChiVariARIBA-reference-sequences)
* [Quick Start](#quick-start)
* [Installation](#installation)
  * [Required dependencies](#required-dependencies)
  * [Using pip3](#using-pip3)
  * [From Source](#from-source)
  * [Docker](#docker)
  * [Singularity](#singularity)
  * [Debian (testing)](#debian-testing)
  * [Ubuntu](#ubuntu)
  * [Dependencies and environment variables](#dependencies-and-environment-variables)
* [Temporary files](#temporary-files)
* [Usage](#usage)
* [Citation](#citation)

## Introduction
ChiVariARIBA is a tool that identifies genes encoding chitin metabolism pathway components by running local assemblies. It utilises the ARIBA tool that was originally developed to identify antibiotic resistance genes

The input is a FASTA file of reference sequences (we have compiled these chitin metabolism pathway genes carefully from the literature) and paired sequencing reads. ChiVariARIBA reports which of the reference sequences were found, plus detailed information on the quality of the assemblies and any variants between the sequencing reads and the reference sequences.

## Preparing ChiVariARIBA reference sequences
add details here of how the earlier work was done - for consistency with the manuscript


We used a quality-controlled collection of 241 annotated genome sequences to create the pangenome used as input for this workflow. Details of how these genomes were collated are available in the manuscript linked to this repository (see below).
From these annotated assemblies, a list of 189 genes associated with chitin metabolism were identified. The list included genes for which experimental or functional evidence was available in the literature supporting their role in chitin metabolism (see linked manuscript). Representative gene sequences from the pangenome gene families which contained these genes were collated into the ariba_chitin_genes.fasta file available in this repository. The fasta file header for each sequence corresponds to the gene family name in the pangenome output, and the annotation and accession number for the representative gene sequence is available in the ariba_metadata.tsv file in this repository. The metadata file also captures all of the sequence variation contained in the sequence alignment for each gene family. When ARIBA is used with this database, the presence and absence of all 189 gene families in read files will be reported, as well as any known sequence variants for genes previously-identified in the original pangenome. 

## Quick Start
We have attached a pre-prepared database to this repository (ariba_database) which can be pulled down and ran independently with your selected input sequencing reads. However, we have also included the original .fasta and .tsv metadata file that were originally used to create ariba_database. You can recreate your own ariba_database using these files with the following command:

Prepare reference data for ChiVariARIBA:

    ariba prepareref -f ariba_chitin_genes.fasta -m ariba_metadata.tsv ariba_database

Run local assemblies and call variants:

    ariba run ariba_database reads1.fastq reads2.fastq ariba_output

Collapse the output report so that each gene found represents one single line (please move the collapse_report.py script into the ariba_output directory before executing):

    python3 collapse_report.py
    
This produces ariba_collapsed_report.tsv which should contain all chitin metabolism pathway genes present in the input reads.

It is also possible to summarise data from several runs:

    ariba summary out.summary ariba_output1/ariba_collapsed_report1.tsv ariba_output2/ariba_collapsed_report2.tsv ariba_output2/ariba_collapsed_report2.tsv

Please read the [ARIBA wiki page][ARIBA wiki] for full usage instructions for ARIBA.

## Tutorials
[The Jupyter notebook tutorial](https://github.com/sanger-pathogens/pathogen-informatics-training)

## Installation

In order to use ChiVariARIBA, you must first install ARIBA. If you encounter an issue when installing ARIBA please contact your local system administrator. If you encounter a bug you can log it [here](https://github.com/sanger-pathogens/ariba/issues).

### Required dependencies
  * [Python3][python] version >= 3.6.0
  * [Bowtie2][bowtie2] version >= 2.1.0
  * [CD-HIT][cdhit] version >= 4.6
  * [MUMmer][mummer] version >= 3.23

ARIBA also depends on several Python packages, all of which are available
via pip. Installing ARIBA with pip3 will get these automatically if they
are not already installed:
  * dendropy >= 4.2.0
  * matplotlib>=3.1.0
  * pyfastaq >= 3.12.0
  * pysam >= 0.9.1
  * pymummer >= 0.10.1
  * biopython

### Using pip3
Install ARIBA using pip:

    pip3 install ariba

### From Source
Download the latest release from this github repository or clone it. Run the tests:

    python3 setup.py test

**Note for OS X:** The tests require gawk which will need to be installed separately, e.g. via Homebrew.

If the tests all pass, install:

    python3 setup.py install

Alternatively, install directly from github using:

    pip3 install git+https://github.com/sanger-pathogens/ariba.git #--user

### Docker
ARIBA can be run in a Docker container. First install Docker, then install the latest
version of ARIBA:

    docker pull gchr.io/sanger-pathogens/ariba:latest

All Docker images are listed in the
[packages page](https://github.com/sanger-pathogens/ariba/pkgs/container/ariba).

To use ARIBA use a command like this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:

    docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/ariba ariba -h

When calling Ariba via Docker (as above) you'll also need to add **/data/** in front of all the passed in file or directory names (e.g. /data/my_output_folder).


### Singularity

ARIBA can be run in a Singularity container. First install Singularity.
[Releases](https://github.com/sanger-pathogens/ariba/releases) include
a Singularity image to download.

Alternatively, build your own Singularity image:

```
singularity build ariba.simg Singularity.def
```


### Debian (Ariba version may not be the latest)
ARIBA is available in the latest version of Debian, and over time will progressively filter through to Ubuntu and other distributions which use Debian. To install it as root:

    sudo apt-get install ariba

### Ubuntu
You can use `apt-get` (see above), or to ensure you get the latest version of ARIBA, the following commands can be
used to install ARIBA and its dependencies. This was tested on a new instance of Ubuntu 16.04.

    sudo  apt-get update
    sudo apt-get install -y python3-dev python3-pip python3-tk zlib1g-dev bowtie2 mummer cd-hit
    export ARIBA_CDHIT=cdhit-est
    sudo pip3 install ariba

### Dependencies and environment variables

By default, ARIBA will look for the dependencies in your `$PATH`, using
the names in the table below. This behaviour can be overridden and
point ARIBA to a specific program using environment variables.
The environment variable is checked first and is used if it is set.
Otherwise ARIBA looks in your `$PATH` for the default name. This applies
to the following dependencies.

| Dependency     |  Default executable    | Environment variable name |
|----------------|------------------------|---------------------------|
| Bowtie2        | `bowtie2`              | `$ARIBA_BOWTIE2`          |
| CD-HIT (est)   | `cd-hit-est`           | `$ARIBA_CDHIT`            |


For example, you could specify an exact version of a bowtie2 executable
that you compiled and downloaded in your home directory (assuming BASH):

    export ARIBA_BOWTIE2=$HOME/bowtie2-2.1.0/bowtie2

Note that ARIBA also runs `bowtie2-build`, for which it uses the
`bowtie2` executable with `-build` appended. So in this case
it would try to use

    $HOME/bowtie2-2.1.0/bowtie2-build

## Temporary files

ARIBA can temporarily make a large number of files whilst running, which
are put in a temporary directory made by ARIBA.  The total size of these
files is small, but there can be a many of them. This can be a
problem when running large numbers (100s or 1000s) of jobs simultaneously
on the same file system.
The parent directory of the temporary directory is determined in the
following order of precedence:

1. The value of the option `--tmp_dir` (if that option was used)
2. The environment variable `$ARIBA_TMPDIR` (if it is set)
3. The environment variable `$TMPDIR` (if it is set)
4. If none of the above is found, then use the run's output directory.

Each temporary directory
is unique to one run of ARIBA, and is automatically deleted at the end
of the run (even if ARIBA was killed by the user or crashed).
For example,

    export $ARIBA_TMPDIR=/tmp

will result in the creation of a new directory inside `/tmp`, which
will have a name of the form

    /tmp/ariba.tmp.abcdef

where the suffix `abcdef` is a random string of characters, chosen
such that `/tmp/ariba.tmp.abcdef` does not already exist.

The exception to the above is if the option `--noclean` is used.
This forces the temporary directory to be placed in the output
directory, and temporary files are kept. It is intended for
debugging.

## Usage
    usage: ariba <command> <options>

    optional arguments:
      -h, --help      show this help message and exit

    Available commands:

	aln2meta      Converts multi-aln fasta and SNPs to metadata
	expandflag    Expands flag column of report file
	flag          Translate the meaning of a flag
	getref        Download reference data
	micplot       Make violin/dot plots using MIC data
	prepareref    Prepare reference data for input to "run"
	pubmlstget    Download species from PubMLST and make db
	pubmlstspecies
		      Get list of available species from PubMLST
	refquery      Get cluster or sequence info from prepareref output
	run           Run the local assembly pipeline
	summary       Summarise multiple reports made by "run"
	test          Run small built-in test dataset
	version       Get versions and exit

Please read the [ARIBA wiki page][ARIBA wiki] for full usage instructions.

## Citation
ChiVariARIBA:A modular, reproducible, and editable workflow for characterising chitin gene variation in Vibrios and related species
Naughton E, Dorman M.J. (Journal)(year)(doi)

ARIBA: rapid antimicrobial resistance genotyping directly from sequencing reads
Hunt M, Mather AE, Sánchez-Busó L, Page AJ, Parkhill J , Keane JA, Harris SR.
Microbial Genomics 2017. doi: [110.1099/mgen.0.000131](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000131)


  [ariba biorxiv]: http://biorxiv.org/content/early/2017/04/07/118000
  [bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  [cdhit]: http://weizhongli-lab.org/cd-hit/
  [ARIBA wiki]: https://github.com/sanger-pathogens/ariba/wiki
  [mummer]: http://mummer.sourceforge.net/
  [python]: https://www.python.org/


