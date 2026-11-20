library(DESeq2)

# sample metadata (like groupfile.txt)
coldata <- read.table("results/leafcutter2_test/groups_file.txt",
                      header = FALSE, col.names = c("sample", "condition"),
                      stringsAsFactors = TRUE)
rownames(coldata) <- coldata$sample

# Col2 (unstranded). Skip 4 summary rows for each sample
star_files <- file.path("results/gene_counts", paste0(coldata$sample, "_ReadsPerGene.out.tab"))

counts_list <- lapply(star_files, function(f) {
  df <- read.table(f, skip = 4, header = FALSE, row.names = 1)
  df[, 1]
})

count_matrix <- do.call(cbind, counts_list)
colnames(count_matrix) <- coldata$sample
rownames(count_matrix) <- rownames(read.table(star_files[1], skip = 4, header = FALSE, row.names = 1))

# Create DESeqDataSet and filter out low-count genes
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

#  differential expression (scz vs bip)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Schizophrenia", "Bipolar"))
res_ordered <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(res_ordered), file = "results/deseq2_scz_vs_bip.csv")
