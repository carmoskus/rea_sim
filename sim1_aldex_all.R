args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

denom = "all"
##test = 1 # 1 for Welch's t-test; 2 for Wilcoxon test

countsN = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

## Convert counts to integer mode manually so I can fix overflows
counts = as.integer(countsN)
counts[is.na(counts)] = .Machine$integer.max
dim(counts) = dim(countsN)
rownames(counts) = rownames(countsN)
colnames(counts) = colnames(countsN)
rm(countsN)

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

#counts = counts[apply(counts, 1, var) > 0,] # Filter out genes with zero variance

## Do DEX
library(ALDEx2)

ald = aldex(counts, as.character(col.info$group), denom=denom)

## Make output data frame
mk.out = function (test) {
    if (test == 1) {
        p = ald$we.ep
    } else {
        p = ald$wi.ep
    }

    df = data.frame(rowMeans(counts[rownames(ald),]), ald$diff.btw, ald$effect, ncol(counts)-2, p)
    colnames(df) = c("baseMean", "log2FC", "stat", "df", "p.value")

    df = df[order(df$p.value),]

    skipped = setdiff(rownames(counts), rownames(ald))
    df2 = data.frame(baseMean=rowMeans(counts[skipped,]))
    df2$log2FC = NA
    df2$stat = NA
    df2$df = NA
    df2$p.value = NA

    df.full = rbind(df, df2)

    name = paste0("aldex_", denom, "_", test)
    write.csv(df.full, file=paste0(subdir, name, "_res.csv"))
}

sapply(1:2, mk.out)

#write.csv(log2(nc+1), file=paste0(subdir, name, "_log2counts.csv"))
#write.table(1/nfs, file=paste0(subdir, name, "_sizes.txt"), sep="\t")
#mc = data.frame(AveLogCPM=y$AveLogCPM, trended.dispersion=y$trended.dispersion, tagwise.dispersion=y$tagwise.dispersion)
#rownames(mc) = rownames(y)
#write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")

