args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: ana_X.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "voom_quantile"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Do voom-limma
library(limma)

design = model.matrix(~group, data=col.info)
v = voom(counts, design, normalize.method="quantile")
fit = lmFit(v, design)
fit = eBayes(fit)

# Make output data frame
df = data.frame(rowMeans(v$E), fit$coefficients[,"groupb"], fit$t[,"groupb"], fit$df.residual, fit$p.value[,"groupb"])
colnames(df) = c("baseMean", "log2FC", "t", "df", "p.value")

df = df[order(df$p.value),]

write.csv(df, file=paste0(subdir, name, "_res.csv"))
