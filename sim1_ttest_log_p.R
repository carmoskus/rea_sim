args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = paste0("ttest_log_", args[3], "_", args[4], "_", args[5])

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, args[3], "_", args[4], ".txt"), header=TRUE, row.names=1, sep="\t")

## Subset the samples
counts = counts[, col.info$subset == args[5]]
col.info = col.info[col.info$subset == args[5],]

## Normalize to total reads
num.reads = colSums(counts)
sfs = num.reads / mean(num.reads)
counts = t(t(counts) / sfs)

## Log-transform counts
counts = log2(counts + 1)

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
