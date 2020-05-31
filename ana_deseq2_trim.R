args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: ana_X.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

countsN = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

## Convert counts to integer mode manually so I can fix overflows
counts = as.integer(countsN)
counts[is.na(counts)] = .Machine$integer.max
dim(counts) = dim(countsN)
rownames(counts) = rownames(countsN)
colnames(counts) = colnames(countsN)
rm(countsN)

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Start DESeq2
library("DESeq2")

se = SummarizedExperiment(counts)
coldata = colData(se)
coldata$Group = col.info$group

dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~ Group)
dds = DESeq(dds)

# Output results
res = results(dds, c('Group', 'b', 'a'))
res = res[order(res$pvalue),]
write.csv(as.data.frame(res), file=paste0(subdir, "deseq2_res.csv"))
