args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]

# Generate row data
n = 1000

p.ex = 1-0.75
m.ex = 9.5
sd.ex = 2.5
m.uex = -0.7
sd.uex = 3.7

ex = rbinom(n, 1, p.ex)
log2means = ifelse(ex, rnorm(n, mean=m.ex, sd=sd.ex), rnorm(n, mean=m.uex, sd=sd.uex))

psi.shape = 2
psi.rate = 6.3
psi.offset = 0.21
psis = rgamma(n, shape=psi.shape, rate=psi.rate)+psi.offset

# Generate counts
ns = 100
x = replicate(ns, rnbinom(n, mu=2^log2means, size=1/psis^2))

# Output
subdir = paste0("a/", arg1, "/")
dir.create(subdir)

# Output created sample
rownames(x) = paste0("G",1:nrow(x))
colnames(x) = paste0("S",1:ncol(x))
write.table(x, file=paste0(subdir, "counts.txt"), quote=FALSE, sep="\t")

# Output meta information
meta = matrix(c("p.ex", p.ex, "m.ex", m.ex, "sd.ex", sd.ex, "m.uex", m.uex, "sd.uex", sd.uex,
                    "psi.shape", psi.shape, "psi.rate", psi.rate, "psi.offset", psi.offset), ncol=2, byrow=TRUE)
write.table(meta, file=paste0(subdir, "meta.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# Output row info
row.info = cbind(log2mean=log2means, psi=psis)
rownames(row.info) = rownames(x)
write.table(row.info, file=paste0(subdir, "rows.txt"), quote=FALSE, sep="\t")

# Output col info
col.info = cbind(group=sample(rep(c("a","b"), ncol(x)/2), ncol(x)))
rownames(col.info) = colnames(x)
write.table(col.info, file=paste0(subdir, "cols.txt"), quote=FALSE, sep="\t")
