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

library(limma)
library(edgeR)

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

name = paste0("lmFit_", arg.norm)

norm = norms[[arg.norm]]
subdir = paste0("sims/", arg.dir, "/")
res.out = paste0("/", name, "_res.csv")

x = lapply(arg.start:arg.end, function (arg.num) {
    counts = as.matrix(read.table(paste0(subdir, arg.num, "/counts.txt"), header=TRUE, sep="\t", row.names=1))
    ## TODO: possibly add following line to all ana2 scripts; if so, also include warning about this
    ## counts[is.na(counts)] = .Machine$integer.max
    col.info = read.table(paste0(subdir, arg.num, "/cols.txt"), header=TRUE, row.names=1, sep="\t")

    dge = DGEList(counts=counts)
    dge = norm(dge)

    eff.lib.sizes = dge$samples$lib.size * dge$samples$norm.factors
    counts =  t(log2(t(counts + 0.5)/(eff.lib.sizes + 1) * 1e+06))

    design = model.matrix(~group, data=col.info)
    
    fit = lmFit(counts, design)
    fit = eBayes(fit)
    
    ## Make output data frame
    df = data.frame(baseMean=rowMeans(counts), log2FC=fit$coefficients[,"groupb"], t=fit$t[,"groupb"], df=fit$df.residual, p.value=fit$p.value[,"groupb"])
    df = df[order(df$p.value),]
    write.csv(df, file=paste0(subdir, arg.num, res.out))
})
