args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]

subdir = paste0("a/", arg1, "/")

name = "edgeR"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Start edgeR
library("edgeR")

mod = model.matrix(~ col.info$group)

y = DGEList(counts=counts)
y = calcNormFactors(y)
y = estimateGLMCommonDisp(y, mod)
y = estimateGLMTrendedDisp(y, mod)
y = estimateGLMTagwiseDisp(y, mod)
fit = glmFit(y, mod)

# Output data
lrt = glmLRT(fit, coef=2)
topTags(lrt, n=50)
write.csv(topTags(lrt, n=100000), file=paste0(subdir, name, "_res.csv"))

#sfs = y$samples$lib.size / mean(y$samples$lib.size)
nfs = y$samples$norm.factors
names(nfs) = colnames(y)
nc = t(t(counts) * nfs)
write.csv(log2(nc+1), file=paste0(subdir, name, "_log2counts.csv"))

# Output metadata
write.table(1/nfs, file=paste0(subdir, name, "_sizes.txt"), sep="\t")

mc = data.frame(AveLogCPM=y$AveLogCPM, trended.dispersion=y$trended.dispersion, tagwise.dispersion=y$tagwise.dispersion)
rownames(mc) = rownames(y)
write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")

