args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]

if (is.na(arg.dir) || nchar(arg.dir) == 0) {
    write("Usage: prog.R subdir analysis", stderr())
    quit(save="no", status=1)
}

analyses = data.frame(name = c("deseq2", "edgeR", "ttest_log", "voom_TMM"),
                      col = c("pvalue", "PValue", "p.value", "p.value"),
                      stringsAsFactors=FALSE)

## Function to check every DEX gene and report whether it is significant at the given threshold, along with its parameters
check = function (i) {
    ## Load row information and pull out DEX genes
    rows = read.table(paste0("sims/", arg.dir, "/", i, "/rows.txt"), header=TRUE, sep="\t", row.names=1)
    rows = rows[rows$log2FC != 0,]

    ## Load analyses
    for (i in 1:nrow(analyses)) {
        analysis = analyses$name[i]
        col = analyses$col[i]

        d = read.csv(paste0("sims/", arg.dir, "/", i, "/", analysis, "_res.csv"), row.names=1)
        d = d[rownames(rows),]

        ## Add p-value to rows
        rows[,analysis] = d[,col]
    }
    rows
}

out = Reduce(rbind, lapply(1:1000, check))
write.table(out, file=paste0("sims/", arg.dir, "/pow2.txt"), sep="\t", quote=FALSE)
