#Import os module for os.system
import os

#Import entrez to access NCBI databases
from Bio import Entrez

#Import SeqIO for parsing files
from Bio import SeqIO

#Import csv for blast output
import csv

#email for Entrez
Entrez.email = "tscott1999@sbcglobal.net"

#Check if the testdata is being used
#usingTestData will be true if so
usingTestData = os.path.isdir('testdata')

#Create the output folder
os.system('mkdir miniProject_Thomas_Scott')

#Copy the required files into the miniProject folder
#(sleuth R script, table for sleuth, and betaherpesvirinae
#fasta file)
os.system('cp miniProject_sleuth.R miniProject_Thomas_Scott')
os.system('cp miniProject_table.txt miniProject_Thomas_Scott')
os.system('cp betaherpes.fasta miniProject_Thomas_Scott')

#if using the test data, copy the data into the folder as well
if usingTestData:
    os.system('cp -r testdata/ miniProject_Thomas_Scott')

#change directory to the miniProject folder
os.chdir('miniProject_Thomas_Scott')

#Create the log file that will
#be used throughout the pipeline
logFile = open('miniProject.log','w')

#make a data directory
os.system('mkdir data')

samples={'Donor 1 (2dpi)':'SRR5660030.1',
         'Donor 1 (6dpi)':'SRR5660033.1',
         'Donor 3 (2dpi)':'SRR5660044.1',
         'Donor 3 (6dpi)':'SRR5660045.1'}

#If the test data is not available, download the full
#data files from SRA and use fastq-dump to uncompress it
if not usingTestData:

    #command and link common to all downloads
    downloadCommand = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/'

    #change directory to the data folder
    os.chdir('data')
    
    #Retrieve the transcriptome data from SRA using wget
    for sample in samples.values():
        os.system(downloadCommand + sample[:-2] + '/' + sample)

    #command for fastq-dump of paired-end reads
    fastqDumpCommand = 'fastq-dump -I --split-files '

    #use fastq dump to uncompress the data and
    #split into paired end read files
    for sample in samples.values():
        os.system(fastqDumpCommand + sample)

    #Exit the data directory
    os.chdir('..')

#if using test data, move the test data to the data directory
else:
    os.system('mv testdata/* data')
    os.system('rm -r testdata/')

#Retrieve the record for EF99921(HCMV) from NCBI (GenBank) using Entrez
handle = Entrez.efetch(db='nucleotide', id='EF999921',rettype='gb',retmode='text')

#Parse the GenBank format using SeqIO.
#records[0] corresponds to the record for EF99921.
records = list(SeqIO.parse(handle,'genbank'))
record = records[0]

#Open a file to write the CDS features to
outfile = open('EF99921_CDS.txt','w')

#Initialize cds to 0 (for the log file)
cds = 0

#This for loop extracts the CDS features from the record,
#writing the protein id, record name (EF99921), and sequence
#to the CDS.txt file. The loop also increments cds for
#each cds feature found.
for feature in record.features:
    if feature.type=='CDS':
        cds += 1 #increment cds for the log file

        #write the protein id, record name (EF99921),
        #and sequence to the CDS.txt file
        outfile.write(">%s %s\n%s\n" % (
            feature.qualifiers['protein_id'][0],
            record.name,
            feature.location.extract(record).seq))

#Close the CDS.txt file
outfile.close()

#Write the entire record to a fasta file for a
#later step in the pipeline (Bowtie2 filtering)
outfile = open('EF99921_Genome.fasta','w')
SeqIO.write(records, outfile, 'fasta')
outfile.close()

#Write the number of CDS features to the log file
logFile.write('The HCMV genome (EF999921) has ' + str(cds)+ ' CDS.\n\n')

#Create a transcriptome index for HCMV using the created
#EF99921_CDS.txt file. Put the index in an index directory.
os.system('mkdir index')
os.system('time kallisto index -i index/index.idx EF99921_CDS.txt')

#Quantify TPM for the CDS features in each transcriptome
#with kallisto. The results will go into a results directory.
os.system('mkdir results')

#kallistoCommand will return the UNIX command for kallisto quant given a sample ID
def kallistoCommand(sampleID):
    return 'time kallisto quant -i index/index.idx -o results/' + sampleID[:-2] + ' -b 30 -t 2 data/' + sampleID + '_1.fastq data/' + sampleID + '_2.fastq'

#Run kallisto quantification for each sample
for sample in samples.values():
    os.system(kallistoCommand(sample))

#Run the sleuth R script to check for differential expression
#between the two timepoints, writing details for significant
#transcripts to the log file.
os.system('Rscript miniProject_sleuth.R')

for line in open('sleuthResults.txt'):
    logFile.write(line)
logFile.write('\n')

#Use Bowtie2 to filter for reads that map to HCMV
#First, create an index for HCMV, using the EF99921_Genome.fasta
#file created earlier, naming the index 'HCMV'
os.system('bowtie2-build EF99921_Genome.fasta HCMV')

#Next, map the reads for each sample to the index.

#bowtie2Command returns the UNIX command used for bowtie2
#mapping given a sampleID, saving the reads that map to the index.
def bowtie2Command(sampleID):
    return 'bowtie2 --quiet -x HCMV -1 data/' + sampleID + '_1.fastq -2 data/' + sampleID + '_2.fastq -S ' + sampleID + 'map.sam --al-conc ' + sampleID + '_mapped_%.fq'

#Run bowtie2 for each sample
for sample in samples.values():
    os.system(bowtie2Command(sample))

#countFilteredReads takes as input a sample's name
#and ID, and counts the number of reads in that sample's
#original fastq file and the mapped fastq file.
#These numbers are then written to the log file.
def countFilteredReads(sampleName,sampleID):

    #Open original data file
    infileBefore=open('data/'+sampleID+'_1.fastq')

    #Open mapped reads file
    infileAfter=open(sampleID+'_mapped_1.fq')

    #Parse both files for the records/reads
    recordsBefore = SeqIO.parse(infileBefore,'fastq')
    recordsAfter = SeqIO.parse(infileAfter,'fastq')

    #Count the number of reads/records in each file
    countBefore = 0
    countAfter = 0
    for record in recordsBefore:
        countBefore+=1
    for record in recordsAfter:
        countAfter+=1

    #Write the counts to the log file
    logFile.write(sampleName + ' had ' + str(countBefore) + ' read pairs before Bowtie2 filtering and '\
                  + str(countAfter) + ' read pairs after.\n')

    #Close the files
    infileBefore.close()
    infileAfter.close()

#Count the reads before and after bowtie2 filtering for each sample
for sampleName, sampleID in samples.items():
    countFilteredReads(sampleName,sampleID)
    
logFile.write('\n')

#Assemble all 4 mapped transcriptomes using SPAdes, and write
#the SPAdes command used to the log file
spadesCommand = 'spades -k 55,77,99,127 -t 2 --only-assembler'\
                ' --pe1-1 SRR5660030.1_mapped_1.fq --pe1-2 SRR5660030.1_mapped_2.fq'\
                ' --pe2-1 SRR5660033.1_mapped_1.fq --pe2-2 SRR5660033.1_mapped_2.fq'\
                ' --pe3-1 SRR5660044.1_mapped_1.fq --pe3-2 SRR5660044.1_mapped_2.fq'\
                ' --pe4-1 SRR5660045.1_mapped_1.fq --pe4-2 SRR5660045.1_mapped_2.fq'\
                ' -o transcriptome_assembly/'
os.system(spadesCommand)
logFile.write('Spades command used for assembly: ' + spadesCommand + '\n\n')

#From the assembly contigs.fasta file, count the number
#of contigs with length>1000

#Open the contigs file, located in the transcriptome_assembly directory
contigsFile = open('transcriptome_assembly/contigs.fasta')

#Parse the fasta file to separate each contig
contigs = SeqIO.parse(contigsFile,'fasta')

#Initialize a count for the number of long contigs (length>1000)
numLongContigs = 0

#Initialze a list for the long contigs
longContigs = []

#For each contig in the file, if its length is > 1000,
#increment the count for the number of long contigs and
#append it to the long contigs list.
for contig in contigs:
    if len(str(contig.seq)) > 1000:
        numLongContigs += 1
        longContigs.append(contig)

#Write the number of long contigs (length>1000) to the logfile
logFile.write('There are ' + str(numLongContigs) + ' contigs > 1000 bp in the assembly.\n')

#From the longContigs, count the number of bp in all
#the contigs.
basePairs = 0
for longContig in longContigs:
    basePairs += len(str(longContig.seq))

#Write the number of bp in contigs>1000 to the log file
logFile.write('There are ' + str(basePairs) + ' bp in the assembly.\n\n')

#Retreive the longest contig (first contig in the longContigs list)
#and write it to its own fasta file.
longestContig = longContigs[0]
SeqIO.write(longestContig, 'longest.fasta', 'fasta')

#Create a local blast+ database of sequences from the
#Betaherpesvirinae subfamily using sequences
#included in the Github repo (from NCBI)
os.system('makeblastdb -in betaherpes.fasta -out betaherpes -title betaherpes -dbtype nucl')

#Use the longest contig from the assembly as input to blast+
#against the local database, creating a blastResults.txt file
#that is tab-delimitted and holds the top hits with the
#appropriate information for each (see headers)
os.system('blastn -query longest.fasta -db betaherpes -out blastResults.txt -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"')

#Write the appropriate headers to the logfile
headers=['sacc','pident','length','qstart','qend','sstart','send','bitscore','evalue','stitle']
for header in headers:
    logFile.write(header + '\t')
logFile.write('\n')

#Open the blast results file and write the top
#10 hits to the logfile
blast_results=open('blastResults.txt')
lines = blast_results.readlines()
for i in range(0,10):
    try:
        logFile.write(lines[i])
    except IndexError:
        break
