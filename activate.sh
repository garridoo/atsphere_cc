#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

### QIIME

export PATH=/projects/dep_psl/grp_psl/garridoo/16S/qiime/bin/:$PATH

export LD_LIBRARY_PATH=/projects/dep_psl/grp_psl/garridoo/16S/qiime/lib/python2.7/site-packages/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/projects/dep_psl/grp_psl/garridoo/16S/qiime/lib/:$LD_LIBRARY_PATH

export PYTHONPATH=/projects/dep_psl/grp_psl/garridoo/16S/qiime/lib/python2.7/site-packages/:$PYTHONPATH

### usearch
export PATH=/projects/dep_psl/grp_psl/garridoo/16S/usearch:$PATH

