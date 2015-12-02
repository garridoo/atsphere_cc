#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# exits whenever a function returns 1
set -e

# get path to scripts
scripts_dir=$(dirname $0)

# load config file
source $scripts_dir/config.sh

# load functions
source $scripts_dir/16s.functions.sh

# activate QIIME, etc.
source $scripts_dir/activate.sh

# parameters

min_length=315
min_id=97

# representative sequences

IPL_rep_seqs=$working_dir/IPL/rep_seqs.fasta
LNS_rep_seqs=$working_dir/leaf_natural_sites/rep_seqs.fasta
sanger_rep_seqs=$data_dir/Leaf_isolates_Sanger.fasta
WGS_rep_seqs=$data_dir/Leaf_isolates_WGS.fasta

# prepare results directory

results_dir=$working_dir/recovery_rates_leaf
mkdir -p $results_dir
rm -f -r $results_dir/*

output=$results_dir/"output.txt"
logfile=$results_dir/"log.txt"

# add label to seq IDs

log "preparing representative OTU sequences..."

cat $IPL_rep_seqs | sed 's/>/>IPL_/g' >> $results_dir/IPL_rep_seqs.fasta
cat $LNS_rep_seqs | sed 's/>/>LNS_/g' >> $results_dir/LNS_rep_seqs.fasta
cat $sanger_rep_seqs | sed 's/>/>sanger_/g' >> $results_dir/sanger_rep_seqs.fasta
cat $WGS_rep_seqs | sed 's/>/>WGS_/g' >> $results_dir/WGS_rep_seqs.fasta

cat $results_dir/IPL_rep_seqs.fasta \
    $results_dir/sanger_rep_seqs.fasta \
    $results_dir/WGS_rep_seqs.fasta \
    >> $results_dir/rec_rep_seqs.fasta

# build database using the representative sequences of the
# leaf-associated OTUs from the culture independent studies

target=$results_dir/rec_rep_seqs.fasta
query=$results_dir/LNS_rep_seqs.fasta

log "building blast database..."

formatdb -n $results_dir/blastDB_leaf \
         -i $target \
         -p false \
         &>> $output

# blast the 16S sequences found in the genomes to
# the reference DB of representative sequences

log "blasting sequences..."

blastall -p blastn \
         -e 1e-10 \
         -m 8 \
         -a $n_cores \
         -v 100000 \
         -b 100000 \
         -i $query \
         -d $results_dir/blastDB_leaf \
         -o $results_dir/blast_results_leaf.txt \
         &>> $output

# cleanup
rm -f blastDB_leaf* formatdb.log error.log

# parse the blast results taking matches > 97% identity
# over at least 99% of the reference sequences length (315 bp)

log "parsing results..."

$scripts_dir/recovery_rates_leaf.R $working_dir \
                                   $results_dir \
                                   $min_length \
                                   $min_id \
                                   &>> $output

log "DONE!"

