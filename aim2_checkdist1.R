args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "aim2_checkdist1"

points = 1:99 / 100
test.ln = function (params) {
    qlnorm(points, meanlog=params[1], sdlog=params[2])
}
test.nb = function (params) {
    qnbinom(points, mu=params[1], size=params[2])
}

d2 = 2^as.matrix(read.csv(paste0(subdir,"edgeR_log2counts.csv"), row.names=1))-1
means = rowMeans(d2)
d2 = d2[means >= quantile(means, p=0.75),]

means = rowMeans(d2)
vars = apply(d2, 1, var)

options(warn = -1)

## Test for log-normal fit
ln.p = data.frame(meanlog = log(means)-1/2*log(vars/means^2+1), sdlog = sqrt(log(vars/means^2+1)))
ln.q = t(apply(ln.p, 1, test.ln))

ks.ln = sapply(1:nrow(d2), function (i) {
    ks.test(ln.q[i,], d2[i,])$p.value
})

## Test for negative binomial fit
nb.p = data.frame(mu=means, size=(means^2)/(vars-means))
nb.q = t(apply(nb.p, 1, test.nb))

ks.nb = sapply(1:nrow(d2), function (i) {
    ks.test(nb.q[i,], d2[i,])$p.value
})

#summary(ks.ln)
#summary(ks.nb)

t = wilcox.test(ks.ln, ks.nb, conf.int=TRUE)
#t
write.table(data.frame(estimate=t$estimate, p.value=t$p.value), file=paste0(subdir, name, ".txt"), row.names=FALSE, quote=FALSE, sep="\t")
