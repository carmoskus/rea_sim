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

name = "ttestH_log"
if (arg.norm != "ms") {
    name = paste0(name, "_", arg.norm)
}

norm = norms[[arg.norm]]
subdir = paste0("sims/", arg.dir, "/")
res.out = paste0("/", name, "_res.csv")

x = lapply(arg.start:arg.end, function (arg.num) {
    counts = as.matrix(read.table(paste0(subdir, arg.num, "/counts.txt"), header=TRUE, sep="\t", row.names=1))
    ## TODO: possibly add following line to all ana2 scripts; if so, also include warning about this
    ## counts[is.na(counts)] = .Machine$integer.max
    col.info = read.table(paste0(subdir, arg.num, "/cols.txt"), header=TRUE, row.names=1, sep="\t")

    ## Normalize and log-transform counts into logCPM with an offset
    dge = DGEList(counts=counts)
    dge = norm(dge)

    eff.lib.sizes = dge$samples$lib.size * dge$samples$norm.factors
    counts =  t(log2(t(counts + 0.5)/(eff.lib.sizes + 1) * 1e+06))
    
    ## Calculate stats for groups 'a' and 'b'
    a.mask = col.info$group=='a'
    b.mask = ! a.mask
    a.n = sum(a.mask)
    b.n = sum(b.mask)
    a = counts[,a.mask]
    b = counts[,b.mask]
    a.m = rowMeans(a)
    b.m = rowMeans(b)
    a.m2 = rowSums(a^2)
    b.m2 = rowSums(b^2)
    a.v = 1/(a.n-1) * (a.m2 - a.n*a.m^2)
    b.v = 1/(b.n-1) * (b.m2 - b.n*b.m^2)

    ## Calculate test stats
    d.m = a.m - b.m
    df = a.n + b.n - 2
    v.p = ((a.n-1)*a.v + (b.n-1)*b.v) / df
    t = d.m / sqrt(v.p*(1/a.n + 1/b.n))
    p = 2*pt(abs(t), df, lower.tail=FALSE)

    ## Make output
    res = data.frame(mean=rowMeans(counts), log2FC=-d.m, t, df, p.value=p)
    res = res[order(res$p.value),]
    
    write.csv(res, file=paste0(subdir, arg.num, res.out))
})
