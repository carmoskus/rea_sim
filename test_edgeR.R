##

arg.dir = "nullA_v5nb2"
arg.num = "1"

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "edgeR"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

## Start edgeR
library(edgeR)

mod = model.matrix(~ group, data=col.info)

y = DGEList(counts=counts)
y = calcNormFactors(y)
y = estimateGLMCommonDisp(y, mod)
y = estimateGLMTrendedDisp(y, mod)
y = estimateGLMTagwiseDisp(y, mod)
fit = glmFit(y, mod)

## DEX gene info
lrt = glmLRT(fit, coef=2)
tt = topTags(lrt, n=100000, sort.by="none")$table
## write.csv(topTags(lrt, n=100000), file=paste0(subdir, name, "_res.csv"))

## edgeR normalization functions
e.cpm = cpm(y)
e.lcpm = cpm(y, log=TRUE)

## Manually normalize counts
nfs = y$samples$norm.factors
eff.lib.sizes = y$samples$lib.size * nfs
sfs = eff.lib.sizes / mean(eff.lib.sizes)
names(sfs) = colnames(y)
nc = t(t(counts) / sfs)
nc.l = log2(nc+1)
## write.csv(log2(nc+1), file=paste0(subdir, name, "_log2counts.csv"))

## Metadata
## write.table(sfs, file=paste0(subdir, name, "_sizes.txt"), sep="\t")

mc = data.frame(AveLogCPM=y$AveLogCPM, trended.dispersion=y$trended.dispersion, tagwise.dispersion=y$tagwise.dispersion)
rownames(mc) = rownames(y)
## write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")

e.m = matrix(c(log2(rowMeans(counts)+1), rowMeans(log2(counts+1)), rowMeans(nc.l), log2(rowMeans(nc)+1), y$AveLogCPM, rowMeans(e.lcpm), log2(rowMeans(e.cpm)+1), tt$logCPM),
             nrow=nrow(counts))
colnames(e.m) = c("logMeanCounts", "meanLogCounts", "nc.l", "log2.nc", "AveLogCPM", "e.lcpm", "log2.cpm", "logCPM")

e.m = e.m[, c("meanLogCounts", "nc.l", "e.lcpm", "logMeanCounts", "log2.nc", "log2.cpm", "AveLogCPM")] ## "logCPM"/tt$logCPM is the exact same as AveLogCPM
e.me = e.m[rowMeans(e.m) > quantile(rowMeans(e.m), probs=0.75),]

e.mec = cor(e.me)

