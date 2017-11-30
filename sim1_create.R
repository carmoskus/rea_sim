args = commandArgs(trailingOnly=TRUE)
arg1 = args[1]

# Generate row data
n = 1000
n.dex = 50
n.zeros = 100
n.tot = n + n.dex + n.zeros

p.ex = 1-0.75
m.ex = 9.5
sd.ex = 2.5
m.uex = -0.7
sd.uex = 3.7

min.fc = 0.05
max.fc = 0.5

ex = c(rbinom(n, 1, p.ex), rep(1, n.dex))
log2means = ifelse(ex, rnorm(n+n.dex, mean=m.ex, sd=sd.ex), rnorm(n+n.dex, mean=m.uex, sd=sd.uex))

psi.shape = 2
psi.rate = 6.3
psi.offset = 0.21
psis = rgamma(n+n.dex, shape=psi.shape, rate=psi.rate)+psi.offset


# Generate counts
#ns = 100
#means = c(2^log2means, rep(0, n.zeros))
#log2means = c(log2means, rep(NA, n.zeros))
#psis = c(psis, rep(1, n.zeros))
#x = replicate(ns, rnbinom(n+n.zeros, mu=means, size=1/psis^2))

# Generate counts for 1 group at a time
ns.g = 50
ns = ns.g * 2
means = c(2^log2means, rep(0, n.zeros))
log2means = c(log2means, rep(NA, n.zeros))
psis = c(psis, rep(1, n.zeros))

# Generate DEXness
group = sample(rep(c("a","b"), ns.g), ns)
log2FC = c(rep(0, n), sample(c(-1,1), n.dex, replace=TRUE)*runif(n.dex, min=min.fc, max=max.fc), rep(0, n.zeros))

# Group a is the base state
a = replicate(ns.g, rnbinom(n.tot, mu=means, size=1/psis^2))

# Group b is altered
b = replicate(ns.g, rnbinom(n.tot, mu=means*2^log2FC, size=1/psis^2))

# Put them back together
x = matrix(0, nrow=n.tot, ncol=ns)
x[, group == "a"] = a
x[, group == "b"] = b

# Output
subdir = paste0("sims/a/", arg1, "/")
dir.create(subdir)

# Output created sample
rownames(x) = paste0("G",1:nrow(x))
colnames(x) = paste0("S",1:ncol(x))
write.table(x, file=paste0(subdir, "counts.txt"), quote=FALSE, sep="\t")

# Output meta information
meta = matrix(c("p.ex", p.ex, "m.ex", m.ex, "sd.ex", sd.ex, "m.uex", m.uex, "sd.uex", sd.uex,
    "psi.shape", psi.shape, "psi.rate", psi.rate, "psi.offset", psi.offset,
    "n", n, "n.zeros", n.zeros, "n.dex", n.dex, "n.tot", n.tot,
    "ns.g", ns.g, "ns", ns,
    "min.fc", min.fc, "max.fc", max.fc), ncol=2, byrow=TRUE)
write.table(meta, file=paste0(subdir, "meta.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# Output row info
row.info = cbind(log2mean=log2means, mean=means, psi=psis, log2FC=log2FC)
rownames(row.info) = rownames(x)
write.table(row.info, file=paste0(subdir, "rows.txt"), quote=FALSE, sep="\t")

# Output col info
col.info = cbind(group=group)
rownames(col.info) = colnames(x)
write.table(col.info, file=paste0(subdir, "cols.txt"), quote=FALSE, sep="\t")
