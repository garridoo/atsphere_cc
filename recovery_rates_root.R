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

SI_otu_table.file <- paste(working_dir, "/start_inoculum/otu_table.txt", sep="")
SI_taxonomy.file <- paste(working_dir, "/start_inoculum/taxonomy.txt", sep="")
blast_results.file <- paste(results_dir, "/blast_results_root.txt", sep="")

# load start inoculum OTU table and taxonomy infomration

SI_taxonomy <- read.table(SI_taxonomy.file, sep="\t", header=F, fill=T)
SI_otu_table <- read.table(SI_otu_table.file, sep="\t", header=T)
rownames(SI_otu_table) <- SI_otu_table[, 1]
SI_otu_table <- SI_otu_table[, -1]

# remove soil samples

idx <- !grepl("InoIRL5and6|InoIRL4and7", colnames(SI_otu_table))
SI_otu_table <- SI_otu_table[, idx]

# remove non-bacterial and Chloroflexi OTUs

SI_taxonomy <- SI_taxonomy[SI_taxonomy[, 2]=="k__Bacteria", ]
SI_taxonomy <- SI_taxonomy[SI_taxonomy[, 3]!="p__Chloroflexi", ]
idx <- rownames(SI_otu_table) %in% SI_taxonomy[, 1]
SI_otu_table <- SI_otu_table[idx, ]

# get top 100 OTUs and those with >.1% relative abundance

SI_otu_table_norm <- apply(SI_otu_table, 2, function(x) x/sum(x)*100)

idx <- apply(SI_otu_table_norm, 1, mean) > .1
#~ idx <- rowSums(SI_otu_table_norm > .1) == dim(SI_otu_table_norm)[2]
#~ idx <- rowSums(SI_otu_table_norm > .1) > 0

abu_otus <- rownames(SI_otu_table)[idx]
abu_otus <- gsub("^", "SI_", abu_otus)

top_100 <- names(sort(apply(SI_otu_table_norm, 1, mean), decreasing=T)[1:100])
top_100 <- gsub("^", "SI_", top_100)

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
cat("root current study:\n")
cat(paste("top 100 OTUs:               ", format(recovery_rate_100, digits=4), "%\n", sep=""))
cat(paste("abundant (> 0.1% R.A) OTUs: ", format(recovery_rate_abu, digits=4), "%\n", sep=""))
cat("----------------------------------\n")
sink()

