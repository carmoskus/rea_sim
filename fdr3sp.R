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
    exp.genes = intersect(exp.genes, rownames(de))
    de = de[exp.genes,]
    nam = names(de)

    if ("p12" %in% nam) {
        exp = de$p12
        dir = de$z12
    } else if ("logCPM" %in% nam) {
        ## edgeR
        exp = de$PValue
        dir = de$logFC
    } else if ("lfcSE" %in% nam) {
        ## DESeq2
        exp = de$pvalue
        dir = de$log2FoldChange
    } else if ("mean" %in% nam) {
        ## t-test
        exp = de$p.value
        dir = de$log2FC
    } else {
        ## Other = voom, aim2
        exp = de$p.value
        dir = de$log2FC
    }

    if (adj != "un") {
        exp = p.adjust(exp, method=adj)
    }
    
    mask = exp <= th
    ns = sum(mask)
    
    if (is.na(ns) || ns == 0) {
        c(NA, 0)
    } else {
        sig = rownames(de)[mask]
        ds = sig[sig %in% dex.genes]

        tde = length(dex.genes)
        
        dir.real = rows[ds, "log2FC"]
        dir.seen = dir[mask][sig %in% dex.genes]

        nds = sum(sign(dir.real) == sign(dir.seen))
        c(1 - nds / ns, nds / tde) # (FDR, power)
    }
}

N = 1000
out = t(sapply(1:N, fdr))
colnames(out) = c("eFDR", "Power")
write.table(out, file=paste0("sims/", arg.dir, "/", analysis, "_fdrp", th, "_", adj, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)

warnings()
