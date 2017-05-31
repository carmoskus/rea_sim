args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]

subdir = paste0("b/", arg1, "/")

name = "edgeR"

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Normalize to total reads
num.reads = colSums(counts)
sfs = num.reads / mean(num.reads)
counts = t(t(counts) / sfs)

# Do t-tests
test = function (row) {
    r = t.test(row ~ col.info$group)
    x = c(mean=mean(row), log2FC=r$estimate[2]-r$estimate[1], t=r$statistic, df=r$parameter, p.value=r$p.value)
    names(x) = c("mean", "log2FC", "t", "df", "p.value")
    x
}
res = as.data.frame(t(apply(counts, 1, test)))

res = res[order(res$p.value),]

# Output
write.csv(res, file=paste0(subdir, "ttest_res.csv"))

