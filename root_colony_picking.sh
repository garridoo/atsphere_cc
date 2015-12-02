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

# prepare results directory

working_dir=$working_dir/root_colony_picking
mkdir -p $working_dir
rm -f -r $working_dir/*

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

### 454 root colony picking (library 866) ###

# parameters

l=866                               # library ID
bc_length=40                        # barcode length
phred=20                            # min. phred score
qual=25                             # min. qual score
bc_err=0                            # allowed barcode errors
t_length=379                        # trimming length

# initialize lib. results directory

rm -f -r $working_dir/"$l"
mkdir -p $working_dir/"$l"

# truncating reads to equal length

log "["$l"] truncating reads..."
truncate_fasta_qual_files.py -f $data_dir/"$l".fasta \
                             -q $data_dir/"$l".qual \
                             -b $t_length \
                             -o $working_dir/"$l"/trunc \
                             &>> $output

mv $working_dir/"$l"/trunc/"$l"_filtered.fasta $working_dir/"$l"/filtered.fasta
mv $working_dir/"$l"/trunc/"$l"_filtered.qual $working_dir/"$l"/filtered.qual
rm -f -r $working_dir/"$l"/trunc

# quality filtering and demultiplexing

log "["$l"] demultiplexing..."
split_libraries.py -f $working_dir/"$l"/filtered.fasta \
                   -q $working_dir/"$l"/filtered.qual \
                   -m $data_dir/"$l"_mapping.txt \
                   -s $qual \
                   -e $bc_err \
                   -b $bc_length \
                   -l $t_length \
                   -o $working_dir/"$l"/demux \
                   &>> $output

mv $working_dir/"$l"/demux/* $working_dir/"$l"
rm -f -r $working_dir/"$l"/demux

# edit barcode label identifier for usearch compatibility
cat $working_dir/"$l"/seqs.fna | \
    sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' \
    >> $working_dir/seqs.fasta

# check sample sizes
sampleSizes $working_dir/seqs.fasta \
            $data_dir/"$l"_mapping.txt \
            $working_dir/sample_sizes.txt \
            &>> $output

# dereplication
log "dereplicating..."
usearch -derep_fulllength $working_dir/seqs.fasta \
        -fastaout $working_dir/seqs_unique.fasta \
        -sizeout \
        &>> $output

# abundance sort and discard singletons
log "sorting by abundance and discarding simgletons..."
usearch -sortbysize $working_dir/seqs_unique.fasta \
        -fastaout $working_dir/seqs_unique_sorted.fasta \
        -minsize $min_size \
        &>> $output

# OTU clustering
log "OTU clustering using UPARSE..."
usearch -cluster_otus $working_dir/seqs_unique_sorted.fasta \
        -otus $working_dir/otus.fasta \
        &>> $output

# chimera detection
log "removing chimeras..."
usearch -uchime_ref $working_dir/otus.fasta \
        -db $gold_db \
        -strand plus \
        -nonchimeras $working_dir/otus_nc.fasta \
        -threads $n_cores \
        &>> $output

# align sequences to database using PyNAST and remove remaining
log "aligning OTU representative sequences to database..."
align_seqs.py -i $working_dir/otus_nc.fasta \
              -t $gg_core_aligned_db \
              -p $min_id_aln \
              -o $working_dir

# rename OTUs and remove alignment gaps

log "renaming OTUs..."

sed -i 's/-//g' $working_dir/otus_nc_aligned.fasta &>> $output

cat $working_dir/otus_nc_aligned.fasta | \
    awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' \
    >> $working_dir/rep_seqs.fasta

# generate OTU table
log "generating OTU table..."
usearch -usearch_global $working_dir/seqs.fasta \
        -db $working_dir/rep_seqs.fasta \
        -strand plus \
        -id $id_threshold \
        -uc $working_dir/read_mapping.uc \
        &>> $output

# convert uc file to txt
log "converting uc OTU table file into text format..."
python $usearch_dir/uc2otutab.py $working_dir/read_mapping.uc \
    1> $working_dir/otu_table.txt \
    2>> $output

# taxonomy assignment
log "taxonomy assignment..."
assign_taxonomy.py -i $working_dir/rep_seqs.fasta \
                   -r $gg_core_db \
                   -t $gg_taxonomy \
                   -m $tax_method \
                   -o $working_dir/tax \
                   &>> $output

# cleanup
mv $working_dir/tax/rep_seqs_tax_assignments.txt $working_dir/taxonomy.txt
rm -f -r $working_dir/tax

sed -i 's/; /\t/g' $working_dir/taxonomy.txt

# convert OTU table to biom
log "converting OTU table to QIIME compatible biom format..."
biom convert -i $working_dir/otu_table.txt \
             -o $working_dir/otu_table.biom \
             --table-type="OTU table" \
             --to-json \
             &>> $output

# align the representative sequences
log "aligning representative sequences..."
clustalo --seqtype=DNA \
         --threads=$n_cores \
         -i $working_dir/rep_seqs.fasta \
         -o $working_dir/rep_seqs_aligned.fasta \
         --full \
         &>> $output

# filter the alignment
filter_alignment.py -i $working_dir/rep_seqs_aligned.fasta \
                    -o $working_dir \
                    &>> $output

# generate tree from alignment using FastTree
log "generating tree..."
make_phylogeny.py -i $working_dir/rep_seqs_aligned_pfiltered.fasta \
                  -o $working_dir/rep_seqs.tree \
                  &>> $output

log "DONE!"
 
