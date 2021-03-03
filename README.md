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

