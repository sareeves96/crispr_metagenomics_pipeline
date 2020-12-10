#!/bin/bash
# main_run

## *********************************
## Dependencies:

# Needs the standalone version of CRISPRCasFinder (CCF).
# Needs the standalone version of blast.
# Makes sure the version of python running is 3.6+.
# Install the Biopython package.
# Make sure the latest version of Perl is running too.

## *********************************


## *********************************
## User inputs

while getopts 'i:c:o:' opt; do
  case "$opt" in
    i) in_file_name=$OPTARG ;;  # name of the input file. 
    c) CCF_dir=$OPTARG ;; # pathway to the CCF directory
    o) out_dir=$OPTARG ;;
  esac
done

path_dir=`pwd`  # the pathway to the current working directory.

CF_dir=$CCF_dir/CasFinder-2.0.3
CCF_pl=$CCF_dir/CRISPRCasFinder.pl
so_file=$CCF_dir/sel392v2.so
out_dir_2=$out_dir/analysis
mkdir $out_dir_2

perl $CCF_pl -in $in_file_name -out $out_dir -cas -cf $CF_dir -keep -def G -html -copyCSS -meta -so $so_file

# extract spacer regions indexed by associated contig_ids, as well as a spacer index
python3 collect_spacers.py --result $out_dir/result.json --out $out_dir_2
# extract cas clusters index by associated contig_id and with each gene
python3 collect_cas.py --result $out_dir/result.json --out $out_dir_2
# blast all input contigs for taxonomy hits
blastn -query $in_file_name -db ref_prok_rep_genomes -outfmt "6 qseqid sseqid scinames scomnames evalue" -evalue 0.001 -max_target_seqs 5 -out $out_dir_2/blast_contigs_vs_prok_rep_genomes.tab -parse_deflines -num_threads 8
# blast all detected putative spacer regions for taxonomy hits
blastn -query $out_dir_2/labelled_spacers.fa -db ref_viruses_rep_genomes -outfmt "6 qseqid sseqid scinames scomnames evalue bitscore pident" -word_size 5 -evalue 0.005 -max_target_seqs 5 -out $out_dir_2/blast_spacers_vs_rep_viruses.tab -parse_deflines -num_threads 8

python3 label_contigs_with_blast.py --contigs $in_file_name --spacers $out_dir_2/labelled_spacers.fa --contigs_blast $out_dir_2/blast_contigs_vs_prok_rep_genomes.tab --spacers_blast $out_dir_2/blast_spacers_vs_rep_viruses.tab --out $out_dir_2 --result $out_dir/result.json

python3 make_cas_class_genus_pie.py --table $out_dir_2/final_output.tab --cas $out_dir_2/Cas_summary.tsv --out $out_dir_2
python3 make_cas_class_genus_pie_alternate.py --table $out_dir_2/final_output.tab --cas $out_dir_2/Cas_summary.tsv --out $out_dir_2
python3 make_genus_species_pie.py --table $out_dir_2/final_output.tab --out $out_dir_2
