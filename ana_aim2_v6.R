args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "aim2_v6"

## Load edgeR and voom-limma with TMM
ed = read.csv(paste0(subdir, "edgeR_res.csv"), row.names=1)
vd = read.csv(paste0(subdir, "voom_TMM_res.csv"), row.names=1)
cd = read.table(paste0(subdir, "aim2_checkdist4.txt"), header=TRUE, sep="\t")

genes = intersect(rownames(ed), rownames(vd))
ed = ed[genes,]
vd = vd[genes,]
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
nd = nd[order(nd$p.value),]

write.csv(nd, paste0(subdir, name, "_res.csv"))
