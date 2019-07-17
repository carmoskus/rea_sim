args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "aim2_v1"

## Function to look at data and return posterior probability of being NB, given it is either NB or log-normal
p.nb = function () {
    0.5
}

## Load edgeR and voom-limma with TMM

## Average test statistics


#write.csv(topTags(lrt, n=100000), file=paste0(subdir, name, "_res.csv"))
