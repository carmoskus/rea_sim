args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.ana = "FC1"

if (is.na(arg.dir) || is.na(arg.ana) || nchar(arg.dir) == 0 || nchar(arg.ana) == 0) {
    write("Usage: prog.R subdir ana", stderr())
    quit(save="no", status=1)
}

checker = function (arg.num) {
    read.table(paste0("sims/", arg.dir, "/", arg.num, "/cor_", arg.ana, ".txt"), header=TRUE, sep="\t")
}

#modes = c("edgeR", "voom_TMM", "ttest_log_TMM", "aim2_v6")

out = Reduce(rbind, lapply(1:1000, checker))

write.table(out, file=paste0("sims/", arg.dir, "/cor_", arg.ana, ".txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

