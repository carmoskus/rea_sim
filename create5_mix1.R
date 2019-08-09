args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

arg.dir = "v5mix1ASHG"
arg.num = "1"

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

root.dir = "sims"
subdir = paste0(root.dir, "/", arg.dir, "/", arg.num, "/")

## FIXME: Check that output files do not already exist
if (file.exists(paste0(subdir, "")) || 
    file.exists(paste0(subdir, "")) || 
    file.exists(paste0(subdir, "")) || 
    file.exists(paste0(subdir, ""))) {
    write("Error: files to be outputted already exist", stderr())
    quit(save="no", status=1)
}

## Read in settings from directory
conf.file = paste0(root.dir, "/", arg.dir, "/meta.txt")
if (!file.exists(conf.file)) {
    write(paste0("Error: no configuration file found at '", conf.file, "'"), stderr())
    quit(save="no", status=1)
}

conf.data = read.table(conf.file, sep="\t", stringsAsFactors=FALSE, row.names=1)
conf.data = rbind(conf.data, mode=6)
conf.data = rbind(conf.data, v=5)
conf = as.list(conf.data$V2)
names(conf) = rownames(conf.data)

## Generate row data
ex = c(rbinom(conf$n, 1, conf$p.ex), rep(1, conf$n.dex))
nm = conf$n+conf$n.dex
n.tot = nm + conf$n.zeros
log2means = ifelse(ex, rnorm(nm, mean=conf$m.ex, sd=conf$sd.ex), rnorm(nm, mean=conf$m.uex, sd=conf$sd.uex))
dist.options = c("create5_lognorm", "create5_nb")
dists = factor(sample(dist.options, n.tot, replace=TRUE), dist.options)

psis = rgamma(nm, shape=conf$psi.shape, rate=conf$psi.rate)+conf$psi.offset

## Generate counts for 1 group at a time
ns = conf$ns.g * 2
means = c(2^log2means, c(0.0000000001, 0)[dists[(nm+1):(n.tot)]])
log2means = c(log2means, rep(NA, conf$n.zeros))
psis = c(psis, rep(1, conf$n.zeros))

## Generate DEXness
group = sample(rep(c("a","b"), conf$ns.g), ns)
sizes = runif(ns, conf$size.min, conf$size.max)
log2FC = c(rep(0, conf$n), sample(c(-1,1), conf$n.dex, replace=TRUE)*runif(conf$n.dex, min=conf$min.fc, max=conf$max.fc), rep(0, conf$n.zeros))

## Generate means incorporating size factors
meansA = matrix(means, nrow=length(means)) %*% sizes[1:conf$ns.g]
meansB = matrix(means*2^log2FC, nrow=length(means)) %*% sizes[conf$ns.g + 1:conf$ns.g]

## Group a is the base state
varsA = meansA + (meansA*psis)^2
meanlogsA = log(meansA)-1/2*log(varsA/meansA^2+1)
varlogsA = log(varsA/meansA^2+1)
## FIXME: rlnorm calls generate warnings because of means=0
a.ln = matrix(rlnorm(conf$ns.g * length(means), meanlog=meanlogsA, sdlog=sqrt(varlogsA)), nrow=length(means))
a.ln = round(a.ln)
a.nb = matrix(rnbinom(conf$ns.g * length(means), mu=meansA, size=1/psis^2), nrow=length(means))

## Group b is altered
varsB = meansB + (meansB*psis)^2
meanlogsB = log(meansB)-1/2*log(varsB/meansB^2+1)
varlogsB = log(varsB/meansB^2+1)
b.ln = matrix(rlnorm(conf$ns.g * length(means), meanlog=meanlogsB, sdlog=sqrt(varlogsB)), nrow=length(means))
b.ln = round(b.ln)
b.nb = matrix(rnbinom(conf$ns.g * length(means), mu=meansB, size=1/psis^2), nrow=length(means))

## Put them back together
x = matrix(0, nrow=n.tot, ncol=ns)
x[as.integer(dists) == 1, group == "a"] = a.ln[as.integer(dists) == 1,]
x[as.integer(dists) == 1, group == "b"] = b.ln[as.integer(dists) == 1,]
x[as.integer(dists) == 2, group == "a"] = a.nb[as.integer(dists) == 2,]
x[as.integer(dists) == 2, group == "b"] = b.nb[as.integer(dists) == 2,]

## Output
dir.create(subdir)

## Output created sample
rownames(x) = paste0("G",1:nrow(x))
colnames(x) = paste0("S",1:ncol(x))
write.table(x, file=paste0(subdir, "counts.txt"), quote=FALSE, sep="\t")

## Output meta information
write.table(conf.data, file=paste0(subdir, "meta.txt"), quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t")

## Output row info
row.info = data.frame(log2mean=log2means, mean=means, psi=psis, log2FC=log2FC, dist=dists)
rownames(row.info) = rownames(x)
write.table(row.info, file=paste0(subdir, "rows.txt"), quote=FALSE, sep="\t")

## Output col info
col.info = data.frame(group=group, stringsAsFactors=FALSE)
col.info[group == "a","size"] = sizes[1:conf$ns.g]
col.info[group == "b","size"] = sizes[conf$ns.g + 1:conf$ns.g]
rownames(col.info) = colnames(x)
write.table(col.info, file=paste0(subdir, "cols.txt"), quote=FALSE, sep="\t")
