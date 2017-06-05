args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]

subdir = paste0("a/", arg1, "/")

name = "pca1"

# Read in data
log2counts = as.matrix(read.csv(paste0(subdir, "deseq2_notrim_log2counts.csv"), row.names=1))

# Filter based on expression
m = rowMeans(log2counts)
sum(m > log2(3.5+1))
log2counts = log2counts[m > log2(3.5+1),]

# Do PCA on data
p = prcomp(t(log2counts))

# Save PCA data
write.csv(p$x, file=paste0(subdir, name, "_values.csv"))

percs = p$sdev^2 / sum(p$sdev^2)
write.table(percs, file=paste0(subdir, name, "_percs.txt"))
