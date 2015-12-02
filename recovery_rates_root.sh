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

IRL_rep_seqs=$working_dir/IRLs/rep_seqs.fasta
RCP_rep_seqs=$working_dir/root_colony_picking/rep_seqs.fasta
SI_rep_seqs=$working_dir/start_inoculum/rep_seqs.fasta
sanger_rep_seqs=$data_dir/Root_isolates_Sanger.fasta
WGS_rep_seqs=$data_dir/Root_isolates_WGS.fasta

# prepare results directory

results_dir=$working_dir/recovery_rates_root
mkdir -p $results_dir
rm -f -r $results_dir/*

output=$results_dir/"output.txt"
logfile=$results_dir/"log.txt"

# add label to seq IDs

log "preparing representative OTU sequences..."

cat $IRL_rep_seqs | sed 's/>/>IRL_/g' >> $results_dir/IRL_rep_seqs.fasta
cat $RCP_rep_seqs | sed 's/>/>RCP_/g' >> $results_dir/RCP_rep_seqs.fasta
cat $SI_rep_seqs | sed 's/>/>SI_/g' >> $results_dir/SI_rep_seqs.fasta
cat $sanger_rep_seqs | sed 's/>/>sanger_/g' >> $results_dir/sanger_rep_seqs.fasta
cat $WGS_rep_seqs | sed 's/>/>WGS_/g' >> $results_dir/WGS_rep_seqs.fasta

# concatenate recovered representative sequences
# from the IRL and colony picking methods

cat $results_dir/IRL_rep_seqs.fasta \
    $results_dir/RCP_rep_seqs.fasta \
    $results_dir/sanger_rep_seqs.fasta \
    $results_dir/WGS_rep_seqs.fasta \
    >> $results_dir/rec_rep_seqs.fasta

# build database using the representative sequence of the
# root-associated OTUs from the culture independent studies

target=$results_dir/rec_rep_seqs.fasta
query=$results_dir/SI_rep_seqs.fasta

log "building blast database..."

formatdb -n $results_dir/blastDB_root \
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
         -d $results_dir/blastDB_root \
         -o $results_dir/blast_results_root.txt \
         &>> $output

# cleanup
rm -f blastDB_root* formatdb.log error.log

# parse the blast results taking matches > 97% identity
# over at least 99% of the reference sequences length (315 bp)

log "parsing results..."

$scripts_dir/recovery_rates_root.R $working_dir \
                                   $results_dir \
                                   $min_length \
                                   $min_id \
                                   &>> $output

log "DONE!"

