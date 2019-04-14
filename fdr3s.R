args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
analysis = args[2]
th = as.numeric(args[3])
adj = args[4]

if (is.na(adj) || nchar(adj) == 0)
    adj = "fdr"

if (is.na(arg.dir) || is.na(th) || nchar(arg.dir) == 0 || is.na(analysis) || nchar(analysis) == 0) {
    write("Usage: prog.R subdir analysis p-threshold [adjustment]", stderr())
    quit(save="no", status=1)
}

## Function to calculate FDR in named run
fdr = function (arg.num) {
    ## Load metadata showing which genes had effects induced
    rows = read.table(paste0("sims/", arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
    exp.genes = rownames(rows)[!is.na(rows$mean) & rows$mean >= 3.5]
    dex.genes = intersect(exp.genes, rownames(rows)[rows$log2FC != 0])
    
    de = read.csv(paste0("sims/", arg.dir, "/", arg.num, "/", analysis, "_res.csv"), row.names=1)
    de = de[exp.genes,]
    nam = names(de)

    if ("logCPM" %in% nam) {
        ## edgeR
        exp = de$PValue
    } else if ("lfcSE" %in% nam) {
        ## DESeq2
        exp = de$pvalue
    } else if ("mean" %in% nam) {
        ## t-test
        exp = de$p.value
    } else {
        ## Other = voom
        exp = de$p.value
    }

    if (adj != "un") {
        exp = p.adjust(exp, method=adj)
    }
    
    sig = rownames(de)[exp <= th]
    ns = length(sig)
    nds = sum(sig %in% dex.genes)
    
    if (ns == 0) {
        NA
    } else {
        1 - nds / ns
    }
}

N = 1000
out = sapply(1:N, fdr)
write(out, paste0("sims/", arg.dir, "/", analysis, "_fdr", th, "_", adj, ".txt"), sep="\t", ncolumns=1)
