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

## Look at how many DEX genes found at FDR < 0.05
## TODO: make it check if output already exists
permute = function (i) {
    function (analysis) {
        ## Correlate the test statistics for the FDR < 0.05 set
        res.d.fn = paste0(subdir, analysis, "_ph2a_", i, "_NA_res.csv")
        res.d = read.csv(res.d.fn, row.names=1)

        nam = names(res.d)
        if ("logCPM" %in% nam) {
            ## edgeR
            p.d = res.d$PValue
            fc.d = res.d$logFC
        } else if ("lfcSE" %in% nam) {
            ## DESeq2
            p.d = res.d$pvalue
            fc.d = res.d$log2FoldChange
        } else if ("mean" %in% nam) {
            ## t-test
            p.d = res.d$p.value
            fc.d = res.d$log2FC
        } else {
            ## Other = voom
            p.d = res.d$p.value
            fc.d = res.d$log2FC
        }

        ## Select FDR < 5% genes in .d
        fdr.d = p.adjust(p.d, method="BH")
        sum(fdr.d < 0.05, na.rm=TRUE)
    }
}

## Look at how many DEX genes found at FDR < 0.05
N = 20
## "deseq2_notrim" "edgeR" "voom_TMM" "voom" "ttest_log_TMM" "ttest_log"
#analyses = c("deseq2_notrim", "edgeR", "voom_TMM", "voom", "ttest_log_TMM", "ttest_log")
analyses = c("edgeR", "voom_TMM", "ttest_log_TMM")

#x = sapply(sapply(analyses, permute), function (f) sapply(1:N, f))
x = t(sapply(sapply(1:N, permute), function (f) sapply(analyses, f)))
write.table(x, file=paste0(subdir, name, "F.txt"), sep="\t", row.names=FALSE, quote=FALSE)
