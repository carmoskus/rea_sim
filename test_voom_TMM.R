##

arg.dir = "nullA_v5nb2"
arg.num = "1"

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "voom_TMM"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Do voom-limma
library(limma)
library(edgeR)

v.dge = DGEList(counts=counts)
v.dge = calcNormFactors(v.dge)

mod = model.matrix(~ group, data=col.info)

v = voom(v.dge, mod)
v.fit = lmFit(v, mod)
v.eb = eBayes(v.fit)

## Make output data frame
v.df = data.frame(rowMeans(v$E), v.eb$coefficients[,"groupb"], v.eb$t[,"groupb"], v.eb$df.residual, v.eb$p.value[,"groupb"])
colnames(v.df) = c("baseMean", "log2FC", "t", "df", "p.value") ## "baseMean" is identical to v.eb$Amean

## df = df[order(df$p.value),]
## write.csv(df, file=paste0(subdir, name, "_res.csv"))

v.m = matrix(c(log2(rowMeans(counts)+1), rowMeans(log2(counts+1)), v.df$baseMean),
             nrow=nrow(counts))
colnames(v.m) = c("logMeanCounts", "meanLogCounts", "voomMean")

v.me = v.m[rowMeans(v.m) > quantile(rowMeans(v.m), probs=0.75),]

v.ec = cor(v.me)

