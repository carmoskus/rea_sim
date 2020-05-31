args = commandArgs(trailingOnly=TRUE)
arg.dir   = args[1]
arg.norm  = args[2]
arg.start = as.integer(args[3])
arg.end   = as.integer(args[4])

if (is.na(arg.dir) || nchar(arg.dir) == 0 || is.na(arg.norm) || nchar(arg.norm) == 0 || is.na(arg.start) || is.na(arg.end) ||
    ! arg.norm %in% c("DE2","TMM","RLE","UQ","ms","ns")) {
    write("Usage: ana2_X.R subdir norm num.start num.end", stderr())
    quit(save="no", status=1)
}

library("edgeR")
library("DESeq2")

norm.TMM = calcNormFactors
norm.RLE = function (x) calcNormFactors(x, method="RLE")
norm.UQ = function (x) calcNormFactors(x, method="upperquartile")
norm.ms = function (x) calcNormFactors(x, method="none")
norm.ns = function (x) {
    y = calcNormFactors(x, method="none")
    y$samples$lib.size = mean(y$samples$lib.size)
    y
}
norms = list(DE2=norm.ms, TMM=norm.TMM, RLE=norm.RLE, UQ=norm.UQ, ms=norm.ms, ns=norm.ns)

name = "d2trim"
if (arg.norm != "DE2") {
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

    dge = DGEList(counts=counts)
    dge = norm(dge)

    eff.lib.sizes = dge$samples$lib.size * dge$samples$norm.factors
    sfs = eff.lib.sizes / mean(eff.lib.sizes)

    ## Start DESeq2
    se = SummarizedExperiment(counts)
    coldata = colData(se)
    coldata$Group = col.info$group

    dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~ Group)
    if (arg.norm != "DE2") {
        sizeFactors(dds) = sfs
    }
    dds = DESeq(dds, minReplicatesForReplace=5)
    
    ## Output results
    res = results(dds, c('Group', 'b', 'a'), cooksCutoff=Inf)
    res = res[order(res$pvalue),]
    write.csv(as.data.frame(res), file=paste0(subdir, arg.num, res.out))
})
