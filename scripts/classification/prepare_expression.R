# Create expression matrix w/ DESeq2
# rm small expression genes
# Normalize for library size and composition
# scale for modelling
# samples x genes matrix for classification
library(DESeq2)

counts <- read.csv(
    "results/classification/expression_counts.csv",
    row.names = 1,
    check.names = FALSE
)

# rm genes w low expression
keep <- colSums(counts >= 10) >= 5
counts <- counts[, keep, drop = FALSE]

# count data
coldata <- data.frame(
    row.names = rownames(counts)
)

dds <- DESeqDataSetFromMatrix(
    countData = round(t(counts)),
    colData = coldata,
    design = ~ 1
)

# stabilize variance
vsd <- varianceStabilizingTransformation(
    dds,
    blind = TRUE
)

# Export samples x genes
expr_matrix <- t(assay(vsd))

write.csv(
    expr_matrix,
    "results/classification/expression_matrix.csv"
)

cat("Expression matrix generated:",
    nrow(expr_matrix), "samples and",
    ncol(expr_matrix), "genes\n")
