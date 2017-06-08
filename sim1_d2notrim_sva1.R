args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]
arg2 = as.integer(args[2])

subdir = paste0("a/", arg1, "/")
name = paste0("d2notrim_sva", arg2)

paste0(subdir, name)

counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Load SVA
sva = read.csv(paste0(subdir, "sva1n", arg2, "_values.csv"), row.names=1)

# Start DESeq2
library("DESeq2")

se = SummarizedExperiment(counts)
coldata = colData(se)
coldata$Group = col.info$group

# Add covariates
form = "Group"
if (!is.na(arg2) & arg2 > 0) {
    for (i in 1:arg2) {
        coldata[,paste0("Covar", i)] = sva[,i]
        form = paste0(form, " + Covar", i)
    }
}

print(form)
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = reformulate(form))
dds = DESeq(dds, minReplicatesForReplace=Inf)

resultsNames(dds)

# Output results
res = results(dds, c('Group', 'a', 'b'), cooksCutoff=Inf)
res = res[order(res$pvalue),]
write.csv(as.data.frame(res), file=paste0(subdir, name, "_res.csv"))

#write.csv(log2(counts(dds, normalized=TRUE)+1), file=paste0(subdir, name, "_log2counts.csv"))

# Output metadata
write.table(sizeFactors(dds), file=paste0(subdir, name, "_sizes.txt"), sep="\t")

mc = as.data.frame(mcols(dds))
rownames(mc) = rownames(dds)
write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")

