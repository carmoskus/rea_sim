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
    y
}
norms = list(TMM=norm.TMM, RLE=norm.RLE, UQ=norm.UQ, ms=norm.ms, ns=norm.ns)

name = "ttest_log"

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
    sfs = eff.lib.sizes / mean(eff.lib.sizes)

    counts = t(t(counts) / sfs)
    
    ## Log-transform counts
    counts = log2(counts+1)
    
    ## Do t-tests
    test = function (row) {
        r = t.test(row ~ col.info$group)
        x = c(mean=mean(row), log2FC=r$estimate[2]-r$estimate[1], t=r$statistic, df=r$parameter, p.value=r$p.value)
        names(x) = c("mean", "log2FC", "t", "df", "p.value")
        x
    }
    res = as.data.frame(t(apply(counts, 1, test)))
    
    res = res[order(res$p.value),]
    
    ## Output
    write.csv(as.data.frame(res), file=paste0(subdir, arg.num, res.out))
})
