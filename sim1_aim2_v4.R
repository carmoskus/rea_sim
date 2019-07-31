args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "aim2_v3"

## Load edgeR and voom-limma with TMM
ed = read.csv(paste0(subdir, "edgeR_res.csv"), row.names=1)
vd = read.csv(paste0(subdir, "voom_TMM_res.csv"), row.names=1)
cd = read.table(paste0(subdir, "aim2_checkdist2.txt"), header=TRUE, sep="\t")

genes = intersect(rownames(ed), rownames(vd))
ed = ed[genes,]
vd = vd[genes,]
cd = cd[genes,]

## Look at GOF data and return posterior probability of being NB, given it is either NB or log-normal
p.nb = ifelse(cd$K == Inf, 1, cd$K/(1+cd$K))

## Average test statistics
new.log2fc = p.nb*ed$logFC + (1-p.nb)*vd$log2FC
new.p = p.nb*ed$PValue + (1-p.nb)*vd$p.value

nd = data.frame(log2FC=new.log2fc, p.value=new.p)
rownames(nd) = genes
nd = nd[order(nd$p.value),]

write.csv(nd, paste0(subdir, name, "_res.csv"))
