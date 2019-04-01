args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.param = args[2]
arg.first = as.numeric(args[3])
## Continues on for as many arguments as the user wants to provide

if (is.na(arg.dir) || is.na(arg.param) || nchar(arg.dir) == 0 || nchar(arg.param) == 0 || is.na(arg.first)) {
    write("Usage: prog.R template param change1 [change2]...", stderr())
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

run = function (p) {
    ## Make an analysis directory with settings as in template, except changing param to p
    cat("Making analysis with", arg.param, "set to", p, "\n")

    p = as.numeric(p)
    conf.data[arg.param,"V2"] = p

    subdir = paste0(root.dir, "/", arg.dir, "_" ,arg.param, "_", p, "/")

    ## Check if dir exists
    if (file.exists(subdir)) {
        cat("ERROR:", subdir, "already exists\n")
        quit(save="no", status=1)
    }
    
    cat("Making", subdir, "\n")
    dir.create(subdir)

    ## Output new meta file
    write.table(conf.data, file=paste0(subdir, "meta.txt"), quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t")

}

out = sapply(args[3:length(args)], run)
