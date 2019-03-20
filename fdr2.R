args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
th = args[2]
adj = args[3]

if (is.na(adj) || nchar(adj) == 0)
    adj = "fdr"

if (is.na(arg.dir) || is.na(th) || nchar(arg.dir) == 0 ||  nchar(th) == 0) {
    write("Usage: prog.R subdir analysis p-threshold", stderr())
    quit(save="no", status=1)
}

## Function to calculate FDR in named run
fdr = function (i) {
    ## Load data
    d = read.table(paste0("sims/", arg.dir, "/", i, "/check2.txt"), header=TRUE, sep="\t", row.names=1)
    function (analysis) {
        if (! analysis %in% names(d)) {
            write(paste0("Could not find analysis '", analysis, "'"), stderr())
            quit(save="no", status=1)
        }

        ## Pull out data we want
        a = d[[analysis]]
        names(a) = rownames(d)
        n = a["n.dex.exp"]
        sg = a[paste0("sde.", adj, th)]
        sa = a[paste0("se.", adj, th)]

        if (sa == 0) {
            return(NA);
        }
        sb = sa - sg
        sb / sa
    }
}

N = 1000
analyses = strsplit("deseq2_notrim   edgeR   voom_TMM        voom    ttest_log_TMM   ttest_log", "\\s+")[[1]]
out = t(sapply(sapply(1:N, fdr), function (f) sapply(analyses, f)))

colnames(out) = sub("^(.*?)\\.se\\..*$", "\\1", colnames(out))

write.table(out, paste0("sims/", arg.dir, "/fdr", th, "_", adj, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
print(summary(out))
##print(t.test(out, mu=as.numeric(th)))
