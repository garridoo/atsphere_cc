#!/opt/share/software/bin/Rscript

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

args <- commandArgs(TRUE)

working_dir <- args[1]
results_dir <- args[2]
min_length <- as.numeric(args[3])
min_id <- as.numeric(args[4])

bulgarelli_otu_table.file <- paste(working_dir, "/bulgarelli/otu_table.txt", sep="")
bulgarelli_taxonomy.file <- paste(working_dir, "/bulgarelli/taxonomy.txt", sep="")
blast_results.file <- paste(results_dir, "/blast_results_bulgarelli.txt", sep="")

# load start inoculum OTU table and taxonomy infomration

bulgarelli_taxonomy <- read.table(bulgarelli_taxonomy.file, sep="\t", header=F, fill=T)
bulgarelli_otu_table <- read.table(bulgarelli_otu_table.file, sep="\t", header=T)
rownames(bulgarelli_otu_table) <- bulgarelli_otu_table[, 1]
bulgarelli_otu_table <- bulgarelli_otu_table[, -1]

# remove soil samples

idx <- !grepl("Bulk", colnames(bulgarelli_otu_table))
bulgarelli_otu_table <- bulgarelli_otu_table[, idx]

# remove non-bacterial and Chloroflexi OTUs

bulgarelli_taxonomy <- bulgarelli_taxonomy[bulgarelli_taxonomy[, 2]=="k__Bacteria", ]
bulgarelli_taxonomy <- bulgarelli_taxonomy[bulgarelli_taxonomy[, 3]!="p__Chloroflexi", ]
idx <- rownames(bulgarelli_otu_table) %in% bulgarelli_taxonomy[, 1]
bulgarelli_otu_table <- bulgarelli_otu_table[idx, ]

# get top 100 OTUs and those with >.1% relative abundance

bulgarelli_otu_table_norm <- apply(bulgarelli_otu_table, 2, function(x) x/sum(x)*100)

idx <- apply(bulgarelli_otu_table_norm, 1, mean) > .1
#~ idx <- rowSums(bulgarelli_otu_table_norm > .1) == dim(bulgarelli_otu_table_norm)[2]
#~ idx <- rowSums(bulgarelli_otu_table_norm > .1) > 0

abu_otus <- rownames(bulgarelli_otu_table)[idx]
abu_otus <- gsub("^", "bulgarelli_", abu_otus)

top_100 <- names(sort(apply(bulgarelli_otu_table_norm, 1, mean), decreasing=T)[1:100])
top_100 <- gsub("^", "bulgarelli_", top_100)

# load blast results

classes <- c("factor", "factor", "numeric", "integer", "integer", "integer",
             "integer", "integer", "integer", "integer", "numeric", "numeric")

blast_results <- read.table(blast_results.file, sep="\t", header=F, colClasses=classes)

blast_results <- blast_results[, c(1, 2, 3, 7, 8)]
colnames(blast_results) <- c("target", "query", "perc_id", "q_start", "q_stop")
blast_results$length <- blast_results$q_stop - blast_results$q_start + 1

# get hits that cover at least min_length of the query sequence at 97% identity

idx <- which(blast_results$perc_id >= min_id & blast_results$length >= min_length)
hits <- unique(blast_results$target[idx])

write.table(gsub(".*_OTU", "OTU", hits),
            paste(results_dir, "/recovered_OTUs.txt", sep=""),
            quote=F, col.names=F, row.names=F)

recovery_rate_abu <- sum(hits %in% abu_otus) * 100 / length(abu_otus)
recovery_rate_100 <- sum(hits %in% top_100) * 100 / length(top_100)

sink(paste(working_dir, "/recovery_rates.txt", sep=""), append=T)
cat("Bulgarelli et al.:\n")
cat(paste("top 100 OTUs:               ", format(recovery_rate_100, digits=4), "%\n", sep=""))
cat(paste("abundant (> 0.1% R.A) OTUs: ", format(recovery_rate_abu, digits=4), "%\n", sep=""))
cat("----------------------------------\n")
sink()

