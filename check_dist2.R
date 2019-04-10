args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir arg.num", stderr())
    quit(save="no", status=1)
}

name="check_dist2"

library(DescTools)

check = function (i) {
    subdir = paste0("sims/", arg.dir, "/", i, "/")
    counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))
    
    means = rowMeans(counts)
    vars = apply(counts, 1, var)
    skews = apply(counts, 1, Skew)
    kurts = apply(counts, 1, Kurt)

    out = data.frame(Mean=means, Var=vars, Skew=skews, Kurt=kurts)
    rownames(out) = rownames(counts)

    write.table(out, file=paste0(subdir, name, ".txt"), sep="\t", quote=FALSE)
}

check(arg.num)
