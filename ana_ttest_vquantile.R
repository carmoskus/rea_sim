args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: ana_X.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "ttest_vquantile"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

## Remove rows that are constant
vars = apply(counts, 1, var)
counts = counts[vars > 0,] # TODO: re-add these rows to the results with NA test statistics

## Normalize to total reads
library(limma)
library(edgeR)

dge = DGEList(counts=counts)

design = model.matrix(~group, data=col.info)

# Use voom to quantile normalize and transform
v = voom(dge, design, normalize.method="quantile")

vars.E = apply(v$E, 1, var)
E = v$E[vars.E > 0,]

## Do t-tests
test = function (row) {
    r = t.test(row ~ col.info$group)
    x = c(mean=mean(row), log2FC=r$estimate[2]-r$estimate[1], t=r$statistic, df=r$parameter, p.value=r$p.value)
    names(x) = c("mean", "log2FC", "t", "df", "p.value")
    x
}
res = as.data.frame(t(apply(E, 1, test)))

res = res[order(res$p.value),]

## Output
write.csv(res, file=paste0(subdir, name, "_res.csv"))
