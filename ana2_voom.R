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
library("limma")

norm.TMM = calcNormFactors
norm.RLE = function (x) calcNormFactors(x, method="RLE")
norm.UQ = function (x) calcNormFactors(x, method="upperquartile")
norm.ms = function (x) calcNormFactors(x, method="none")
norm.ns = function (x) {
    y = calcNormFactors(x, method="none")
    y$samples$lib.size = mean(y$samples$lib.size)
    y
}
norms = list(TMM=norm.TMM, RLE=norm.RLE, UQ=norm.UQ, ms=norm.ms, ns=norm.ns)

name = "voom"
if (arg.norm != "ms") {
    name = paste0(name, "_", arg.norm)
}

norm = norms[[arg.norm]]
subdir = paste0("sims/", arg.dir, "/")
res.out = paste0("/", name, "_res.csv")

x = lapply(arg.start:arg.end, function (arg.num) {
    counts = as.matrix(read.table(paste0(subdir, arg.num, "/counts.txt"), header=TRUE, sep="\t", row.names=1))
    col.info = read.table(paste0(subdir, arg.num, "/cols.txt"), header=TRUE, row.names=1, sep="\t")

    dge = DGEList(counts=counts)
    dge = norm(dge)

    mod = model.matrix(~ group, data=col.info)

    v = voom(dge, mod)
    fit = lmFit(v, mod)
    fit = eBayes(fit)

    ## Output data
    df = data.frame(rowMeans(v$E), log2(rowMeans(2^v$E)), fit$coefficients[,"groupb"], fit$t[,"groupb"], fit$df.residual, fit$p.value[,"groupb"])
    colnames(df) = c("baseMean", "logOfMeans", "log2FC", "t", "df", "p.value")
    
    df = df[order(df$p.value),]
    
    write.csv(df, file=paste0(subdir, arg.num, res.out))
})
