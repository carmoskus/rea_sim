args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

root.dir = "sims"
subdir = paste0(root.dir, "/", arg.dir, "/", arg.num, "/")

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

## Function to do one permutation of cross-validation and return
## TODO: make it check if output already exists
permute = function (i) {
    function (analysis) {
        ## Correlate the test statistics for the FDR < 0.05 set
        res.d.fn = paste0(subdir, analysis, "_ph1a_", i, "_d_res.csv")
        res.d = read.csv(res.d.fn, row.names=1)
        res.r.fn = paste0(subdir, analysis, "_ph1a_", i, "_r_res.csv")
        res.r = read.csv(res.r.fn, row.names=1)
        res.r = res.r[rownames(res.d),]

        nam = names(res.r)
        if ("logCPM" %in% nam) {
            ## edgeR
            p.d = res.d$PValue
            p.r = res.r$PValue
            fc.d = res.d$logFC
            fc.r = res.r$logFC
        } else if ("lfcSE" %in% nam) {
            ## DESeq2
            p.d = res.d$pvalue
            p.r = res.r$pvalue
            fc.d = res.d$log2FoldChange
            fc.r = res.r$log2FoldChange
        } else if ("mean" %in% nam) {
            ## t-test
            p.d = res.d$p.value
            p.r = res.r$p.value
            fc.d = res.d$log2FC
            fc.r = res.r$log2FC
        } else {
            ## Other = voom
            p.d = res.d$p.value
            p.r = res.r$p.value
            fc.d = res.d$log2FC
            fc.r = res.r$log2FC
        }

        ## Select FDR < 5% genes in .d
        fdr.d = p.adjust(p.d, method="BH")
        mask = !is.na(fdr.d) & fdr.d < 0.05
        p.d = p.d[mask]
        p.r = p.r[mask]
        fc.d = fc.d[mask]
        fc.r = fc.r[mask]

        ## Convert to z-statistics and return the correlation
        z.d = sign(fc.d) * -qnorm(p.d/2)
        z.r = sign(fc.r) * -qnorm(p.r/2)
        cor(z.d, z.r)
    }
}

## Run the cross-validation metric on N permutations
N = 20
## "deseq2_notrim" "edgeR" "voom_TMM" "voom" "ttest_log_TMM" "ttest_log"
analyses = c("deseq2_notrim", "edgeR", "voom_TMM", "voom", "ttest_log_TMM", "ttest_log")

#x = sapply(sapply(analyses, permute), function (f) sapply(1:N, f))
x = t(sapply(sapply(1:N, permute), function (f) sapply(analyses, f)))
write.table(x, file=paste0(subdir, name, "F.txt"), sep="\t", row.names=FALSE, quote=FALSE)
