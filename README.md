# DomainSeq

## Overview
This project is a pipeline for analysis of long-read enrichment processes. Specifically, if a library is generated composed of multiple domains, each member of each library in each selection cannot be observed easily through Illumina sequencing due to read-length restrictions. Additionally, PacBio provides the necesssary read-length, but is throughput limited. In these cases, multiple sources of sequencing data can be integrated to allow rapid, deep sequencing of libraries. Three types of files need to be provided:

   1. An Excel document linking each possible domain to a DNA sequence
   2. A PacBio .fastq file, containing long-reads connecting a particular barcode to the read of the DNA containing the domain set
   3. An arbitrary number of Illumina sequencing files, containing short-reads which cover barcodes observed in specific stages during selection
   
**DomainSeq** is a pipeline that pieces these three sets of information together and presents processed visuals and documents describing the library state throughout selection. This takes the form:

   Illumina data (BC) -> PacBio data (BC to DNA sequence) -> Excel document (DNA sequence to Domains)
   
A number of user controls are provided for each step of analysis along this analysis and can be scaled to any dataset size, taking appropriate steps to minimize computation time.

### TODO:
  - Multithreading

## Requirements
This package has only been tested on Unix/Linux (OSX, Ubuntu, etc). The project is not platform specific, but installation instructions only apply to these systems. Additionally, all code is formated for Python 3.X.

## Installation
We recommend using this package through a [virtual environment](https://docs.python.org/3/library/venv.html). This creates an "instance" of Python that is isolated from the rest of your system. To do this, open terminal and do the following:

1. Clone/download this package, unzip it, and move into the base directory (".../DomainSeq/") using [change directory](https://www.rapidtables.com/code/linux/cd.html) commands.

2. Install the virtual environment package to your native Python:

   > python3 -m pip install --user virtualenv

3. Use this package to make a virtual environment of Python 3:

   > python3 -m virtualenv domainseq-venv
   
   Note: This creates a folder in your current directory called "domainseq-venv" which contains your virtual environment. Call it whatever you want!
   
4. Activate your virtual environment by declaring it as your source for Python:

   > source domainseq-venv/bin/activate

   Note: This command causes the virtual environment (your "domainseq-venv" folder) to be prioritized as the version of Python your session will use.
   
5. Install the required packages for the project:

   > pip install -r requirements.txt
   
Great! If these commands worked without errors, you have:
+ Copied the project
+ Created a virtual Python environment
+ Install the required packages for the project

See the wiki for information on usage.
