To install python dependencies:

`pip3 install -r requirements.txt`

Requires standalone versions of:

CRISPRCasFinder
https://github.com/dcouvin/CRISPRCasFinder

NCBI Blast+ 2.11.0
`cd $HOME`\n
`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz`\n
`tar -zxvf ncbi-blast-2.11.0+-x64-linux.tar.gz`

For local database searching, requires two BLAST databases:
WARNING: the ref_prok_rep_genomes database is >6GB in size
`mkdir blastdb`
`export BLASTDB=$HOME/blastdb`
`cd blastdb`
`../ncbi-blast-2.11.0+/bin/update_blastdb.pl --decompress ref_prok_rep_genomes`
`../ncbi-blast-2.11.0+/bin/update_blastdb.pl --decompress ref_viruses_rep_genomes`

To run the tool, Blast and CRISPRCasFinder must be properly installed.
The main script, main_run.sh, has 3 arguments:
`sh main_run.sh -c [PATH_TO_CRISPRCasFinder_FOLDER] -i [CONTIG_FASTA_FILE] -o [OUTPUT_FOR_CCF]`

A sample contig fasta file is included in the repository: sample_contigs.fa. The contigs were assembled from paired end reads of a human gut metagenome using megahit. The paired end reads can be found here: 
https://www.ebi.ac.uk/ena/browser/view/SRR341725

The results from the CRISPRCasFinder tool are sent to the -o directory, and inside that directory is an analysis subfolder containing additional analysis conducted by the authors of this pipeline.
