args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.param = args[2]
arg.end = as.numeric(args[3])
arg.steps = as.integer(args[4])

if (is.na(arg.dir) || is.na(arg.param) || nchar(arg.dir) == 0 || nchar(arg.param) == 0 || is.na(arg.end) || is.na(arg.steps)) {
    write("Usage: prog.R template param end steps", stderr())
    quit(save="no", status=1)
}

#subdir = paste0("sims/", arg.dir, "/", arg.num, "/")
#name = "deseq2_notrim"

