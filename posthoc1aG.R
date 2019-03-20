args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]

if (is.na(arg.dir) || nchar(arg.dir) == 0) {
    write("Usage: prog.R subdir", stderr())
    quit(save="no", status=1)
}

root.dir = "sims"
subdir = paste0(root.dir, "/", arg.dir, "/")

name = "ph1a"

conf.file = paste0(root.dir, "/", arg.dir, "/meta.txt")
if (!file.exists(conf.file)) {
    write(paste0("Error: no configuration file found at '", conf.file, "'"), stderr())
    quit(save="no", status=1)
}

conf.data = read.table(conf.file, sep="\t", stringsAsFactors=FALSE, row.names=1)
conf = as.list(conf.data$V2)
names(conf) = rownames(conf.data)

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")
a.names = rownames(col.info)[col.info$group == "a"]
b.names = rownames(col.info)[col.info$group == "b"]

## The number of simulations with permutations run on them
N = 100
## "deseq2_notrim" "edgeR" "voom_TMM" "voom" "ttest_log_TMM" "ttest_log"
analyses = c("deseq2_notrim", "edgeR", "voom_TMM", "voom", "ttest_log_TMM", "ttest_log")

summarize = function (arg.num) {
    function (analysis) {
        
    }
}

x = t(sapply(sapply(1:N, summarize), function (f) sapply(analyses, f)))
write.table(x, file=paste0(subdir, name, "G.txt"), sep="\t", row.names=FALSE, quote=FALSE)
