args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

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
