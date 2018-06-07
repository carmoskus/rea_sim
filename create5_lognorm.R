args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

root.dir = "sims"
subdir = paste0(root.dir, "/", arg.dir, "/", arg.num, "/")

## Check that output files do not already exist
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
conf.data = rbind(conf.data, mode=3)
conf.data = rbind(conf.data, v=4)
conf = as.list(conf.data$V2)
names(conf) = rownames(conf.data)

## Generate row data
ex = c(rbinom(conf$n, 1, conf$p.ex), rep(1, conf$n.dex))
nm = conf$n+conf$n.dex
n.tot = nm + conf$n.zeros
log2means = ifelse(ex, rnorm(nm, mean=conf$m.ex, sd=conf$sd.ex), rnorm(nm, mean=conf$m.uex, sd=conf$sd.uex))

psis = rgamma(nm, shape=conf$psi.shape, rate=conf$psi.rate)+conf$psi.offset

## Generate counts for 1 group at a time
ns = conf$ns.g * 2
means = c(2^log2means, rep(0.0000000001, conf$n.zeros))
log2means = c(log2means, rep(NA, conf$n.zeros))
psis = c(psis, rep(1, conf$n.zeros))

## Generate DEXness
group = sample(rep(c("a","b"), conf$ns.g), ns)
sizes = runif(ns, conf$size.min, conf$size.max)
log2FC = c(rep(0, conf$n), sample(c(-1,1), conf$n.dex, replace=TRUE)*runif(conf$n.dex, min=conf$min.fc, max=conf$max.fc), rep(0, conf$n.zeros))

## Group a is the base state
vars = means + (means*psis)^2
meanlogs = log(means)-1/2*log(vars/means^2+1)
varlogs = log(vars/means^2+1)
a = replicate(conf$ns.g, rlnorm(n.tot, meanlog=meanlogs, sdlog=sqrt(varlogs)))

## Group b is altered
meansb = means*2^log2FC
varsb = meansb + (meansb*psis)^2
meanlogsb = log(meansb)-1/2*log(varsb/meansb^2+1)
varlogsb = log(varsb/meansb^2+1)
b = replicate(conf$ns.g, rlnorm(n.tot, meanlog=meanlogsb, sdlog=sqrt(varlogsb)))

## Alter group values by size factors
a = t(round(t(a)*sizes[1:conf$ns.g]))
a[a < 0] = 0
b = t(round(t(b)*sizes[conf$ns.g + 1:conf$ns.g]))
b[b < 0] = 0

## Put them back together
x = matrix(0, nrow=n.tot, ncol=ns)
x[, group == "a"] = a
x[, group == "b"] = b

## Output
dir.create(subdir)

## Output created sample
rownames(x) = paste0("G",1:nrow(x))
colnames(x) = paste0("S",1:ncol(x))
write.table(x, file=paste0(subdir, "counts.txt"), quote=FALSE, sep="\t")

## Output meta information
write.table(conf.data, file=paste0(subdir, "meta.txt"), quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t")

## Output row info
row.info = cbind(log2mean=log2means, mean=means, psi=psis, log2FC=log2FC)
rownames(row.info) = rownames(x)
write.table(row.info, file=paste0(subdir, "rows.txt"), quote=FALSE, sep="\t")

## Output col info
col.info = data.frame(group=group, stringsAsFactors=FALSE)
col.info[group == "a","size"] = sizes[1:conf$ns.g]
col.info[group == "b","size"] = sizes[conf$ns.g + 1:conf$ns.g]
rownames(col.info) = colnames(x)
write.table(col.info, file=paste0(subdir, "cols.txt"), quote=FALSE, sep="\t")
