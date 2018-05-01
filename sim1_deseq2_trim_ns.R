args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Start DESeq2
library("DESeq2")

se = SummarizedExperiment(counts)
coldata = colData(se)
coldata$Group = col.info$group

dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~ Group)
sizeFactors(dds) = rep(1, ncol(counts))
dds = DESeq(dds)

# Output results
res = results(dds, c('Group', 'a', 'b'))
res = res[order(res$pvalue),]
write.csv(as.data.frame(res), file=paste0(subdir, "deseq2_res.csv"))

write.csv(log2(counts(dds, normalized=TRUE)+1), file=paste0(subdir,"deseq2_log2counts.csv"))

# Output metadata
write.table(sizeFactors(dds), file=paste0(subdir, "deseq2_sizes.txt"), sep="\t")

mc = as.data.frame(mcols(dds))
rownames(mc) = rownames(dds)
write.table(mc, file=paste0(subdir, "deseq2_meta.txt"), sep="\t")

