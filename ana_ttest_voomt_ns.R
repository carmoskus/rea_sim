args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: ana_X.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "ttest_voomt_ns"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

## Normalize by TMM
#num.reads = colSums(counts)
#sfs = num.reads / mean(num.reads)
library(limma)
library(edgeR)

dge = DGEList(counts=counts)
dge = calcNormFactors(dge)

eff.lib.sizes = rep(mean(dge$samples$lib.size), length(dge$samples$lib.size))

## Voom-transform counts
counts =  t(log2(t(counts + 0.5)/(eff.lib.sizes + 1) * 1e+06))

## Check for extremely low variance
vars = apply(counts, 1, var)
counts = counts[vars > 0,]

## Do t-tests
test = function (row) {
    r = t.test(row ~ col.info$group)
    x = c(mean=mean(row), log2FC=r$estimate[2]-r$estimate[1], t=r$statistic, df=r$parameter, p.value=r$p.value)
    names(x) = c("mean", "log2FC", "t", "df", "p.value")
    x
}
res = as.data.frame(t(apply(counts, 1, test)))

res = res[order(res$p.value),]

## Output
write.csv(res, file=paste0(subdir, name, "_res.csv"))
