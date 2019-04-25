args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "deseq2_notrim_RLE"

countsN = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

## Convert counts to integer mode manually so I can fix overflows
counts = as.integer(countsN)
counts[is.na(counts)] = .Machine$integer.max
dim(counts) = dim(countsN)
rownames(counts) = rownames(countsN)
colnames(counts) = colnames(countsN)
rm(countsN)

library(limma)
library(edgeR)

dge = DGEList(counts=counts)
dge = calcNormFactors(dge, method="RLE")

eff.lib.sizes = dge$samples$lib.size * dge$samples$norm.factors
sfs = eff.lib.sizes / mean(eff.lib.sizes)

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Start DESeq2
library("DESeq2")

se = SummarizedExperiment(counts)
coldata = colData(se)
coldata$Group = col.info$group

dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~ Group)
sizeFactors(dds) = sfs
dds = DESeq(dds, minReplicatesForReplace=Inf)

# Output results
res = results(dds, c('Group', 'b', 'a'), cooksCutoff=Inf)
res = res[order(res$pvalue),]
write.csv(as.data.frame(res), file=paste0(subdir, name, "_res.csv"))

write.csv(log2(counts(dds, normalized=TRUE)+1), file=paste0(subdir, name, "_log2counts.csv"))

# Output metadata
write.table(sizeFactors(dds), file=paste0(subdir, name, "_sizes.txt"), sep="\t")

mc = as.data.frame(mcols(dds))
rownames(mc) = rownames(dds)
write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")

