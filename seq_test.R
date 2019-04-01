args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.param = args[2]
arg.end = as.numeric(args[3])
arg.steps = as.integer(args[4])

if (is.na(arg.dir) || is.na(arg.param) || nchar(arg.dir) == 0 || nchar(arg.param) == 0 || is.na(arg.end) || is.na(arg.steps)) {
    write("Usage: prog.R template param end steps", stderr())
    quit(save="no", status=1)
}

root.dir = "sims"

##subdir = paste0("sims/", arg.dir, "/", arg.num, "/")
##name = "deseq2_notrim"

conf.file = paste0(root.dir, "/", arg.dir, "/meta.txt")
if (!file.exists(conf.file)) {
    write(paste0("Error: no configuration file found at '", conf.file, "'"), stderr())
    quit(save="no", status=1)
}

conf.data = read.table(conf.file, sep="\t", stringsAsFactors=FALSE, row.names=1)
conf = as.list(conf.data$V2)
names(conf) = rownames(conf.data)

conf.param = conf[[arg.param]]
cat(arg.steps, "instances for", arg.param, "from", conf.param, "to", arg.end, "\n")

for (p in seq(conf.param, arg.end, length.out=arg.steps)) {
    cat(p, "\n")
}
