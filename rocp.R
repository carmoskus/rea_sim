args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
analysis = args[2]
adj = args[3]

if (is.na(adj) || nchar(adj) == 0)
    adj = "un"

if (is.na(arg.dir) || is.na(analysis) || nchar(arg.dir) == 0 || nchar(analysis) == 0) {
    write("Usage: prog.R subdir analysis adjustment", stderr())
    quit(save="no", status=1)
}


## Function to calculate FDR in named run
roc.exp = function (i) {
    ## Load data
    d = read.table(paste0("sims/", arg.dir, "/", i, "/check4.txt"), header=TRUE, sep="\t", row.names=1)
    if (! analysis %in% names(d)) {
        write(paste0("Could not find analysis '", analysis, "'"), stderr())
        quit(save="no", status=1)
    }

    ## Pull out data we want
    a = d[[analysis]]
    names(a) = rownames(d)
    ne = a["n.all.exp"]
    nde = a["n.dex.exp"]
    nce = ne - nde

    thresholds = grep(paste0("^se.", adj), names(a), value=TRUE)
    se = a[thresholds]
    sde = a[sub("^se", "sde", thresholds)]
    sdae = a[sub("^se", "sdae", thresholds)]
    sce = se - sdae
    fpr = sce / (nce + (sde - sdae))
    tpr = sdae / nde
    fdr = 1 - sdae / (sce + sdae)

    data.frame(FPR=fpr, TPR=tpr, FDR=fdr)
}

out = Reduce(rbind, lapply(1:1000, roc.exp))
colnames(out) = c("FPR", "TPR", "FDR")
write.table(out, file=paste0("sims/", arg.dir, "/", analysis, "_roc_", adj, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
