args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]
arg2 = as.integer(args[2])

subdir = paste0("a/", arg1, "/")
name = paste0("sva1n", arg2)

paste0(subdir, name)

# Read in data
log2counts = as.matrix(read.csv(paste0(subdir, "deseq2_notrim_log2counts.csv"), row.names=1))

# Read in sample data
col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

# Filter based on expression
m = rowMeans(log2counts)
sum(m > log2(3.5+1))
log2counts = log2counts[m > log2(3.5+1),]

# Do SVA on data with the indicated number of SVs
library(sva)
mod = model.matrix(~ group, data=col.info)
mod0 = model.matrix(~1, data=col.info)
svobj = sva(as.matrix(log2counts), mod, mod0, n.sv = arg2)

# Save PCA data
write.csv(svobj$sv, file=paste0(subdir, name, "_values.csv"))

