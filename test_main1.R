##

arg.dir = "nullA_v5nb2"
arg.num = "1"

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "main1"

## Load libraries
library(limma)
library(edgeR)

## Load basic data
counts = as.matrix(read.table(paste0(subdir, "counts.txt"), header=TRUE, sep="\t", row.names=1))

col.info = read.table(paste0(subdir, "cols.txt"), header=TRUE, row.names=1, sep="\t")
mod = model.matrix(~ group, data=col.info)

dge = DGEList(counts=counts)
dge = calcNormFactors(dge)

## Manually normalize counts
eff.lib.sizes = dge$samples$lib.size * dge$samples$norm.factors
sfs = eff.lib.sizes / mean(eff.lib.sizes)
names(sfs) = colnames(dge)
nc = t(t(counts) / sfs)
nc.l = log2(nc+1)
## write.csv(log2(nc+1), file=paste0(subdir, name, "_log2counts.csv"))
## write.table(sfs, file=paste0(subdir, name, "_sizes.txt"), sep="\t"

## -------------------------------------------------------
## Do edgeR stuff
e.y = estimateGLMCommonDisp(dge, mod)
e.y = estimateGLMTrendedDisp(e.y, mod)
e.y = estimateGLMTagwiseDisp(e.y, mod)
e.fit = glmFit(e.y, mod)

## DEX gene info
e.lrt = glmLRT(e.fit, coef=2)
e.tt = topTags(e.lrt, n=100000, sort.by="none")$table
## write.csv(topTags(lrt, n=100000), file=paste0(subdir, name, "_res.csv"))

## edgeR normalization functions
e.cpm = cpm(e.y)
e.lcpm = cpm(e.y, log=TRUE)

e.m = matrix(c(log2(rowMeans(counts)+1), rowMeans(log2(counts+1)), rowMeans(nc.l), log2(rowMeans(nc)+1), rowMeans(e.lcpm), log2(rowMeans(e.cpm)+1), e.tt$logCPM),
             nrow=nrow(counts))
colnames(e.m) = c("logMeanCounts", "meanLogCounts", "nc.l", "log2.nc", "e.lcpm", "log2.cpm", "logCPM")

e.m = e.m[, c("meanLogCounts", "nc.l", "e.lcpm", "logMeanCounts", "log2.nc", "log2.cpm", "logCPM")]
e.me = e.m[rowMeans(e.m) > quantile(rowMeans(e.m), probs=0.75),]

e.mec = cor(e.me)

## --------------------------------------------------------------
## Do voom-limma stuff
v = voom(dge, mod)
v.fit = lmFit(v, mod)
v.eb = eBayes(v.fit)

## Make output data frame
v.df = data.frame(rowMeans(v$E), v.eb$coefficients[,"groupb"], v.eb$t[,"groupb"], v.eb$df.residual, v.eb$p.value[,"groupb"])
colnames(v.df) = c("baseMean", "log2FC", "t", "df", "p.value") ## "baseMean" is identical to v.eb$Amean

## df = df[order(df$p.value),]
## write.csv(df, file=paste0(subdir, name, "_res.csv"))

v.m = matrix(c(log2(rowMeans(counts)+1), rowMeans(log2(counts+1)), v.df$baseMean),
             nrow=nrow(counts))
colnames(v.m) = c("logMeanCounts", "meanLogCounts", "voomMean")

v.me = v.m[rowMeans(v.m) > quantile(rowMeans(v.m), probs=0.75),]

v.ec = cor(v.me)

## --------------------------------------------------------------
## Do BMA stuff
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
        ks.test(nc[i,], plnorm, meanlog=params$meanlog[i], sdlog=params$sdlog[i])$p.value
    }
}
test.nb = function (i) {
    if (params$means[i] == 0) {
        1
    } else {
        ks.test(nc[i,], pmynb, mu=params$means[i], size=params$size[i])$p.value
    }
}

## GOF v4
params = data.frame(means = rowMeans(nc), vars = apply(nc, 1, var))
params$meanlog = log(params$means)-1/2*log(params$vars/params$means^2+1)
params$sdlog = sqrt(log(params$vars/params$means^2+1))
params$size = (params$means^2)/ifelse(params$vars-params$means > 0, params$vars-params$means, 0.1)

options(warn = -1)

## Test for log-normal fit
ks.ln = sapply(1:nrow(nc), test.ln)
ks.ln[ks.ln == 0] = min(ks.ln[ks.ln != 0]) * 1e-3

## Test for negative binomial fit
ks.nb = sapply(1:nrow(nc), test.nb)
ks.nb[ks.nb == 0] = min(ks.nb[ks.nb != 0]) * 1e-3

options(warn = 0)

cd = data.frame(ln.p=ks.ln, nb.p=ks.nb, K=ks.nb/ks.ln)
rownames(cd) = rownames(nc)

## Combine datasets and produce average results for genes that have data in both sets
genes = intersect(rownames(e.tt), rownames(v.df))
ed = e.tt[genes,]
vd = v.df[genes,]
cd = cd[genes,]

## Look at GOF data and return posterior probability of being NB, given it is either NB or log-normal
p.nb = ifelse(cd$K == Inf, 1, cd$K/(1+cd$K))

## Average log2fc
new.log2fc = p.nb*ed$logFC + (1-p.nb)*vd$log2FC

## Estimate standard error by back converting the p-value to a z-value and then multiplying by
ed.z = qnorm(ed$PValue/2, lower.tail=FALSE)
ed.se = abs(ed$logFC) / ed.z

vd.z = qnorm(vd$p.value/2, lower.tail=FALSE)
vd.se = abs(vd$log2FC) / vd.z

new.var = (ed.se^2 + ed$logFC^2)*p.nb + (vd.se^2 + vd$log2FC^2)*(1-p.nb) - new.log2fc^2

new.z = new.log2fc / sqrt(new.var)
new.p = pnorm(abs(new.z), lower.tail=FALSE)*2

nd = data.frame(log2FC=new.log2fc, var=new.var, p.value=new.p)
rownames(nd) = genes

## nd = nd[order(nd$p.value),]
## write.csv(nd, paste0(subdir, name, "_res.csv"))

## --------------------------------------------------------------------
## Look at comparison of all estimates for expression level
