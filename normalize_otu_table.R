#!/opt/share/software/bin/Rscript

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

library(metagenomeSeq, quietly=T, warn.conflicts=F)
library(biom, quietly=T, warn.conflicts=F)

args <- commandArgs(TRUE)

input.file <- args[1]
output.file <- args[2]

b <- load_biom(input.file)
p <- cumNormStatFast(b)
b <- cumNorm(b, p=p)
write_biom(MRexperiment2biom(b, norm=TRUE, log=TRUE), output.file)

