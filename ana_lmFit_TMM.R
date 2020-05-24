args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: ana_X.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "lmFit_TMM"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

## Normalize by TMM
library(limma)
library(edgeR)

dge = DGEList(counts=counts)
dge = calcNormFactors(dge)

eff.lib.sizes = dge$samples$lib.size * dge$samples$norm.factors
sfs = eff.lib.sizes / mean(eff.lib.sizes)

counts = t(t(counts) / sfs)

## Log-transform counts
counts = log2(counts + 1)
## Voom-transform counts
#counts =  t(log2(t(counts + 0.5)/(eff.lib.sizes + 1) * 1e+06))

design = model.matrix(~group, data=col.info)

fit = lmFit(counts, design)
fit = eBayes(fit)

# Make output data frame
df = data.frame(rowMeans(counts), fit$coefficients[,"groupb"], fit$t[,"groupb"], fit$df.residual, fit$p.value[,"groupb"])
colnames(df) = c("baseMean", "log2FC", "t", "df", "p.value")

df = df[order(df$p.value),]

write.csv(df, file=paste0(subdir, name, "_res.csv"))
