Welcome to Thomas Scott's COMP 483 Mini Project!

This project consists of a Python wrapper that automates several software tools used in bioinformatics.
It was made in and is intended to run in a UNIX terminal.
The tools used include kallisto, sleuth, Bowtie2, SPAdes, and BLAST+.
You will need to have these tools installed for the pipeline to run correctly.

A full list of tools to be installed and links to instructions are included here:

kallisto: https://pachterlab.github.io/kallisto/source

R: https://www.r-project.org/

sleuth (R package): https://pachterlab.github.io/sleuth/download

dplyr (R package): https://www.r-project.org/nosvn/pandoc/dplyr.html

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2

SPAdes: https://cab.spbu.ru/files/release3.15.1/manual.html

BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Biopython: https://biopython.org/wiki/Download


This pipeline uses transcriptome data of HCMV (Human cytomegalovirus/Human herpesvirus 5) 2- and 6-days post-infection.
The pipeline produces a variety of outputs and writes important information to a log file titled 'miniProject.log'.

Aside from the Python wrapper ('MiniProject.py'), some other needed files are included.
miniProject_sleuth is the sleuth script ran during the pipeline, with miniProject_table.txt being used in that script.
Furthermore, a betaherpes.fasta file is included, so that a local BLAST+ database can be created for betaherpesvirinae towards the end of the pipeline.

As the pipeline requires a large amount of time to run, test data that is a small subset of reads from the input data is included.
If you would like to use the test data, simply include the testdata folder/directory in the same folder as the python wrapper ('MiniProject.py').
If you would like to run the full pipeline, you may move or delete the testdata folder.

Below are instructions on running the code once the required tools are installed:

1.)
