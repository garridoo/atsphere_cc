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

# process data from the different studies

#~ log "processing IRL data..."
#~ $scripts_dir/IRLs.sh
#~ 
#~ log "processing root colony picking data..."
#~ $scripts_dir/root_colony_picking.sh
#~ 
#~ log "processing IPL data..."
#~ $scripts_dir/IPL.sh
#~ 
#~ log "processing root starting inoculum data..."
#~ $scripts_dir/start_inoculum.sh
#~ 
#~ log "processing Bulgarelli et al. data..."
#~ $scripts_dir/bulgarelli.sh
#~ 
#~ log "processing Schlaeppi et al. data..."
#~ $scripts_dir/schlaeppi.sh
#~ 
#~ log "processing leaf natural sites data..."
#~ $scripts_dir/leaf_natural_sites.sh
#~  
#~ rm -f $working_dir/recovery_rates.txt

# calculate recovery rates

log "calculating root recovery rates..."
$scripts_dir/recovery_rates_root.sh

log "calculating Bulgarelli et al. recovery rates..."
$scripts_dir/recovery_rates_bulgarelli.sh

log "calculating Schlaeppi et al. recovery rates..."
$scripts_dir/recovery_rates_schlaeppi.sh

log "calculating leaf recovery rates..."
$scripts_dir/recovery_rates_leaf.sh

log "DONE!"

