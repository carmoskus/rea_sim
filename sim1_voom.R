args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]

subdir = paste0("b/", arg1, "/")

name = "voom"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Do voom-limma
library(limma)

design = model.matrix(~group, data=col.info)
v = voom(counts, design)
fit = lmFit(v, design)
fit = eBayes(fit)

# Make output data frame
df = data.frame(rowMeans(v$E), fit$coefficients[,"groupb"], fit$t[,"groupb"], fit$df.residual, fit$p.value[,"groupb"])
colnames(df) = c("baseMean", "log2FC", "t", "df", "p.value")

df = df[order(df$p.value),]

write.csv(df, file=paste0(subdir, name, "_res.csv"))

#write.csv(log2(nc+1), file=paste0(subdir, name, "_log2counts.csv"))
#write.table(1/nfs, file=paste0(subdir, name, "_sizes.txt"), sep="\t")
#mc = data.frame(AveLogCPM=y$AveLogCPM, trended.dispersion=y$trended.dispersion, tagwise.dispersion=y$tagwise.dispersion)
#rownames(mc) = rownames(y)
#write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")

