args = commandArgs(trailingOnly=TRUE)
arg.dir   = args[1]
arg.start = as.integer(args[2])
arg.end   = as.integer(args[3])

if (is.na(arg.dir) || nchar(arg.dir) == 0 || is.na(arg.start) || is.na(arg.end)) {
    write("Usage: mix_X.R subdir num.start num.end", stderr())
    quit(save="no", status=1)
}

name = "mix_test1"

path = paste0("sims/", arg.dir, "/")

x = lapply(arg.start:arg.end, function (arg.num) {
    spath = paste0(path, arg.num, "/")
    e.res = read.csv(paste0(spath, "edgeR_TMM_res.csv"), row.names=1)
    v.res = read.csv(paste0(spath, "voom_TMM_res.csv"), row.names=1)
    genes = intersect(rownames(e.res), rownames(v.res))
    e.res[genes,]
    v.res[genes,]
    
    fc1 = e.res$logFC
    fc2 = v.res$log2FC
    z1 = sign(fc1)*qnorm(e.res$PValue/2, lower.tail=FALSE)
    p1 = pnorm(abs(z1), lower.tail=FALSE)*2
    z2 = sign(fc2)*qnorm(v.res$p.value/2, lower.tail=FALSE)
    p2 = pnorm(abs(z2), lower.tail=FALSE)*2
    sd1 = fc1 / z1
    sd2 = fc2 / z2
    sd3 = fc2 / v.res$t

    z12 = fc1 / sd2
    p12 = pnorm(abs(z12), lower.tail=FALSE)*2
    z21 = fc2 / sd1
    p21 = pnorm(abs(z21), lower.tail=FALSE)*2

    z13 = fc1 / sd3
    p13 = pnorm(abs(z13), lower.tail=FALSE)*2

    d = data.frame(z1, p1, z2, p2, z12, p12, z13, p13)
    rownames(d) = genes
    write.csv(d, file=paste0(spath, name, "_res.csv"), quote=FALSE)
})
