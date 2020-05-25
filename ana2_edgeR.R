args = commandArgs(trailingOnly=TRUE)
arg.dir   = args[1]
arg.norm  = args[2]
arg.start = as.integer(args[3])
arg.end   = as.integer(args[4])

if (is.na(arg.dir) || nchar(arg.dir) == 0 || is.na(arg.norm) || nchar(arg.norm) == 0 || is.na(arg.start) || is.na(arg.end) ||
    ! arg.norm %in% c("TMM","RLE","UQ","ms","ns")) {
    write("Usage: ana2_X.R subdir norm num.start num.end", stderr())
    quit(save="no", status=1)
}

library("edgeR")

norm.TMM = calcNormFactors
norm.RLE = function (x) calcNormFactors(x, method="RLE")
norm.UQ = function (x) calcNormFactors(x, method="upperquartile")
norm.ms = function (x) calcNormFactors(x, method="none")
norm.ns = function (x) {
    y = calcNormFactors(x, method="none")
    y$samples$lib.size = mean(y$samples$lib.size)
}
norms = list(TMM=norm.TMM, RLE=norm.RLE, UQ=norm.UQ, ms=norm.ms, ns=norm.ns)

if (arg.norm == "TMM") {
    name = "edgeR"
} else {
    name = paste0("edgeR_", arg.norm)
}

norm = norms[[arg.norm]]

x = lapply(arg.start:arg.end, function (arg.num) {
    subdir = paste0("sims/", arg.dir, "/", arg.num, "/")
    counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))
    col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

    mod = model.matrix(~ col.info$group)

    y = DGEList(counts=counts)
    y = norm(y)
    y = estimateGLMCommonDisp(y, mod)
    y = estimateGLMTrendedDisp(y, mod)
    y = estimateGLMTagwiseDisp(y, mod)
    fit = glmFit(y, mod)
    
    ## Output data
    lrt = glmLRT(fit, coef=2)
    write.csv(topTags(lrt, n=100000), file=paste0(subdir, name, "_res.csv"))

    ## nfs = y$samples$norm.factors
    ## eff.lib.sizes = y$samples$lib.size * nfs
    ## sfs = eff.lib.sizes / mean(eff.lib.sizes)
    ## names(sfs) = colnames(y)
    ## nc = t(t(counts) / sfs)
    ## write.csv(log2(nc+1), file=paste0(subdir, name, "_log2counts.csv"))

    ## Output metadata
    ## write.table(sfs, file=paste0(subdir, name, "_sizes.txt"), sep="\t")

    ## mc = data.frame(AveLogCPM=y$AveLogCPM, trended.dispersion=y$trended.dispersion, tagwise.dispersion=y$tagwise.dispersion)
    ## rownames(mc) = rownames(y)
    ## write.table(mc, file=paste0(subdir, name, "_meta.txt"), sep="\t")
})

