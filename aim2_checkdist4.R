args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "aim2_checkdist4"

## GOF v4
pmynb = function (q, mu, size) {
    q.u = floor(q + 1)
    d = 1 - (q.u - q) # 1 - distance from q to q.u
    p.old = pnbinom(q, mu=mu, size=size)
    p.next = pnbinom(q.u, mu=mu, size=size)
    p.new = p.old + d * (p.next - p.old)
    p.new
}
test.ln = function (i) {
    if (params$means[i] == 0) {
        1
    } else {
        ks.test(d2[i,], plnorm, meanlog=params$meanlog[i], sdlog=params$sdlog[i])$p.value
    }
}
test.nb = function (i) {
    if (params$means[i] == 0) {
        1
    } else {
        ks.test(d2[i,], pmynb, mu=params$means[i], size=params$size[i])$p.value
    }
}

d2 = 2^as.matrix(read.csv(paste0(subdir,"edgeR_log2counts.csv"), row.names=1))-1

params = data.frame(means = rowMeans(d2), vars = apply(d2, 1, var))
params$meanlog = log(params$means)-1/2*log(params$vars/params$means^2+1)
params$sdlog = sqrt(log(params$vars/params$means^2+1))
params$size = (params$means^2)/ifelse(params$vars-params$means > 0, params$vars-params$means, 0.1)

options(warn = -1)

## Test for log-normal fit
ks.ln = sapply(1:nrow(d2), test.ln)
ks.ln[ks.ln == 0] = min(ks.ln[ks.ln != 0]) * 1e-3

## Test for negative binomial fit
ks.nb = sapply(1:nrow(d2), test.nb)
ks.nb[ks.nb == 0] = min(ks.nb[ks.nb != 0]) * 1e-3

out = data.frame(ln.p=ks.ln, nb.p=ks.nb, K=ks.nb/ks.ln)
rownames(out) = rownames(d2)

write.table(out, file=paste0(subdir, name, ".txt"), row.names=TRUE, quote=FALSE, sep="\t")
