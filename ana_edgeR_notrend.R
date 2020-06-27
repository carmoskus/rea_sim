args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: ana_X.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "edgeR_notrend"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Start edgeR
library(edgeR)

mod = model.matrix(~ col.info$group)

y = DGEList(counts=counts)
y = calcNormFactors(y)
y = estimateGLMCommonDisp(y, mod)
y = estimateGLMTagwiseDisp(y, mod)
fit = glmFit(y, mod)

# Output data
lrt = glmLRT(fit, coef=2)
write.csv(topTags(lrt, n=100000), file=paste0(subdir, name, "_res.csv"))
