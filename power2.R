args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
analysis = args[2]
th = args[3]

if (is.na(arg.dir) || is.na(analysis) || is.na(th) || nchar(arg.dir) == 0 || nchar(analysis) == 0 || nchar(th) == 0) {
    write("Usage: prog.R subdir analysis p-threshold", stderr())
    quit(save="no", status=1)
}

## Function to check every DEX gene and report whether it is significant at the given threshold, along with its parameters
check = function (i) {
    ## Load data
}

out = sapply(1:1000, check)
write(out, paste0("sims/", arg.dir, "/", analysis, "_pow2_", th, ".txt"), sep="\t", ncolumns=1)
