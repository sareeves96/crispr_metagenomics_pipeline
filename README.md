![alt text](https://github.com/sareeves96/crispr_metagenomics_pipeline/blob/main/genus_species_pie.png?raw=true)

To install python dependencies:

`pip3 install -r requirements.txt`
***

Requires standalone versions of:

<b> CRISPRCasFinder </b> <br/>
https://github.com/dcouvin/CRISPRCasFinder

<b> NCBI Blast+ 2.11.0 </b> <br/>
`cd $HOME`  
`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz`  
`tar -zxvf ncbi-blast-2.11.0+-x64-linux.tar.gz`

For local database searching, requires two BLAST databases: <br/>
<u> WARNING: the ref_prok_rep_genomes database is >6GB in size </u> 
`mkdir blastdb`  
`export BLASTDB=$HOME/blastdb`  
`cd blastdb`  
`../ncbi-blast-2.11.0+/bin/update_blastdb.pl --decompress ref_prok_rep_genomes`  
`../ncbi-blast-2.11.0+/bin/update_blastdb.pl --decompress ref_viruses_rep_genomes`  

***
To run the tool, Blast and CRISPRCasFinder must be properly installed.  
The main script, main_run.sh, has 3 arguments:  
`sh main_run.sh -c [PATH_TO_CRISPRCasFinder_FOLDER] -i [CONTIG_FASTA_FILE] -o [OUTPUT_FOR_CCF]`
***
A sample contig fasta file is included in the repository: sample_contigs.fa. The contigs were assembled from paired end reads of a human gut metagenome using megahit. The paired end reads can be found here:  
https://www.ebi.ac.uk/ena/browser/view/SRR341725
***
The results from the CRISPRCasFinder tool are sent to the -o directory, and inside that directory is an analysis subfolder containing additional analysis conducted by the authors of this pipeline.
***
<b> Final output files:  </b> <br/>
<br/>
<u> TAB seperated files: </u> <br/>

final_output.tab -> All contigs labelled using blast, and any corresponding putative CRISPR spacer sequences. <br/>
casgene_occurence.tab -> all cas genes detected and their occurences as counts. <br/>
genera_occurence.tab -> all genera identified using blast and their occurences as counts. <br/>
genera_cas_occurence.tab -> all genera containing cas genes and the occurence of genes with respect to each genus. <br/>
genera_cas_fraction.tab -> the fraction of each genera containing cas genes as a percentage. </br>
<br/>

<u> PDF files: </u> <br/>

casgene_distribution.pdf -> a pie_chart showing the distribution of cas genes and a table showing the data from the corresponing .tab file <br/>
genera_stats.pdf -> all genera related analysis. <br/>

***
***
