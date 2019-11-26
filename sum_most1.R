##
## Look at summary statistics pulled from the results files of already performed analyses
## Especially compare different methods for estimating expression and fold change
##

arg.dir = "nullA_v5nb2"
arg.num = "1"

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "sum_most1"

## Load libraries
library(limma)
library(edgeR)

## Load basic data
counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")

e.res = read.csv(paste0(subdir, "edgeR_res.csv"), header=TRUE, row.names=1)
v.res = read.csv(paste0(subdir, "voom_TMM_res.csv"), header=TRUE, row.names=1)
b.res = read.csv(paste0(subdir, "aim2_v6_res.csv"), header=TRUE, row.names=1)
d.res = read.csv(paste0(subdir, "deseq2_notrim_res.csv"), header=TRUE, row.names=1)
t.res = read.csv(paste0(subdir, "ttest_log_TMM_res.csv"), header=TRUE, row.names=1)

## Find subset of genes in all analyses and re-order data
genes = intersect(rownames(e.res), rownames(v.res))
genes = intersect(genes, rownames(b.res))
genes = intersect(genes, rownames(d.res))
genes = intersect(genes, rownames(t.res))
n.s = ncol(counts)

counts = counts[genes,]
logMeanCounts = log2(rowMeans(counts)+1)
exp.mask = logMeanCounts > quantile(logMeanCounts, p=0.75)
genes = genes[exp.mask]
logMeanCounts = logMeanCounts[exp.mask]
n.g = length(genes)

counts = counts[genes,]
counts.l = log2(counts+1)
e.res = e.res[genes,]
v.res = v.res[genes,]
b.res = b.res[genes,]
d.res = d.res[genes,]
t.res = t.res[genes,]

## Derive SD of estimates for log2FC from data
e.res$z = qnorm(e.res$PValue/2, lower.tail=FALSE)
e.res$sd = abs(e.res$logFC) / e.res$z

v.res$z = qnorm(v.res$p.value/2, lower.tail=FALSE)
v.res$sd = abs(v.res$log2FC) / v.res$z

d.res$z = qnorm(d.res$pvalue/2, lower.tail=FALSE)
d.res$sd = abs(d.res$log2FoldChange) / d.res$z

t.res$z = qnorm(t.res$p.value/2, lower.tail=FALSE)
t.res$sd = abs(t.res$log2FC) / t.res$z

## Manually normalize counts
dge = DGEList(counts=counts)
dge = calcNormFactors(dge)

eff.lib.sizes = dge$samples$lib.size * dge$samples$norm.factors
sfs = eff.lib.sizes / mean(eff.lib.sizes)
names(sfs) = colnames(dge)
nc = t(t(counts) / sfs)
nc.l = log2(nc+1)

## Do manual estimates
u.lmfc = apply(counts, 1, function (r) {
    x = tapply(r, col.info$group, mean)
    log2(x[2] / x[1])
})
u.mlfc = apply(counts.l, 1, function (r) {
    x = tapply(r, col.info$group, mean)
    x[2] - x[1]
})
n.lmfc = apply(nc, 1, function (r) {
    x = tapply(r, col.info$group, mean)
    log2(x[2] / x[1])
})
n.mlfc = apply(nc.l, 1, function (r) {
    x = tapply(r, col.info$group, mean)
    x[2] - x[1]
})

u.lms = apply(counts, 1, function (r) {
    x = tapply(r, col.info$group, sd)
    log2(mean(x)+1)
}) / sqrt(n.s)
u.mls = apply(counts.l, 1, function (r) {
    x = tapply(r, col.info$group, sd)
    mean(x)
}) / sqrt(n.s)
n.lms = apply(nc, 1, function (r) {
    x = tapply(r, col.info$group, sd)
    log2(mean(x)+1)
}) / sqrt(n.s)
n.mls = apply(nc.l, 1, function (r) {
    x = tapply(r, col.info$group, sd)
    mean(x)
}) / sqrt(n.s)

## Combine expression estimates
m.e = matrix(c(logMeanCounts, rowMeans(log2(counts+1)), log2(rowMeans(nc)+1), rowMeans(nc.l), e.res$logCPM, v.res$baseMean, v.res$logOfMeans, log2(d.res$baseMean+1), t.res$mean),
             nrow=n.g)
colnames(m.e) = c("logMeanCounts", "meanLogCounts", "logMeanNorm", "meanLogNorm", "e.logCPM", "v.baseMean", "v.logOfMeans", "d.logMean", "t.mean")

m.e = m.e[,c("logMeanCounts", "logMeanNorm", "e.logCPM", "d.logMean", "v.logOfMeans", "meanLogCounts", "meanLogNorm", "v.baseMean", "t.mean")]

## Combine FC estimates
m.f = matrix(c(u.lmfc, u.mlfc, n.lmfc, n.mlfc, e.res$logFC, v.res$log2FC, b.res$log2FC, d.res$log2FoldChange, t.res$log2FC),
             nrow=n.g)
colnames(m.f) = c("u.lmfc", "u.mlfc", "n.lmfc", "n.mlfc", "e.logFC", "v.log2FC", "b.log2FC", "d.log2FC", "t.log2FC")

m.f = m.f[,c("u.lmfc", "n.lmfc", "e.logFC", "d.log2FC", "b.log2FC", "u.mlfc", "n.mlfc", "v.log2FC", "t.log2FC")]

## Combine variance estimates
m.s = matrix(c(u.lms, u.mls, n.lms, n.mls, e.res$sd, v.res$sd, sqrt(b.res$var), d.res$sd, t.res$sd),
             nrow=n.g)
colnames(m.s) = c("u.lms", "u.mls", "n.lms", "n.mls", "e.sd", "v.sd", "b.sd", "d.sd", "t.sd")

## Filter for minimum expression
## m.em = rowMeans(m.e)
## m.mask = m.em > quantile(m.em, probs=0.75)
## m.e = m.e[m.mask,]
## m.f = m.f[m.mask,]
## m.s = m.s[m.mask,]

## Make correlation matrices
c.e = cor(m.e)
c.f = cor(m.f)
c.s = cor(m.s)

c.e
c.f
c.s

## cor(m.e, method="spearman")
## cor(m.f, method="spearman")
## cor(m.s, method="spearman")

