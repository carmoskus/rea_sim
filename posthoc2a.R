args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

root.dir = "sims"
subdir = paste0(root.dir, "/", arg.dir, "/", arg.num, "/")

name = "ph2a"

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

## Function to do one permutation of cross-validation and return
## TODO: make it check if output already exists
permute = function (i) {
        ## Split into two groups, seperately for 'a' and 'b' groups
        p = 0.8
        n1 = p * conf$ns.g
        a.d = sample(a.names, n1)
        a.r = setdiff(a.names, a.d)

        n2 = p * conf$ns.g
        b.d = sample(b.names, n2)
        b.r = setdiff(b.names, b.d)

        ## Save the groups
        d = c(a.d, b.d)
        df = data.frame(group = col.info$group, subset = ifelse(rownames(col.info) %in% d, "d", "r"))
        rownames(df) = rownames(col.info)
        filename = paste0(subdir, name, "_", i, ".txt") 
        write.table(df, file=filename, quote=FALSE, sep="\t")
        
    function (analysis) {
        ## Run analysis on discovery samples
        #cat(paste0("echo Rscript sim1_", analysis, "_p.R ", arg.dir, " ", arg.num, " ", name, " ", i, " d | qsub -cwd\n"))
        cat(paste0("Rscript sim1_", analysis, "_p.R ", arg.dir, " ", arg.num, " ", name, " ", i, " d\n"))
        
        ## Run analysis on replication samples
        #cat(paste0("echo Rscript sim1_", analysis, "_p.R ", arg.dir, " ", arg.num, " ", name, " ", i, " r | qsub -cwd\n"))
        cat(paste0("Rscript sim1_", analysis, "_p.R ", arg.dir, " ", arg.num, " ", name, " ", i, " r\n"))

        ## Correlate the test statistics for the FDR < 0.05 set
        ## TODO: need to do this after the previous commands have finished
    }
}

## Run the cross-validation metric on N permutations
N = 20
## "deseq2_notrim" "edgeR" "voom_TMM" "voom" "ttest_log_TMM" "ttest_log"
#analyses = c("deseq2_notrim", "edgeR", "voom_TMM", "voom", "ttest_log_TMM", "ttest_log")
analyses = c("edgeR", "voom_TMM", "ttest_log_TMM")

#x = sapply(sapply(analyses, permute), function (f) sapply(1:N, f))
x = sapply(sapply(1:N, permute), function (f) sapply(analyses, f))
