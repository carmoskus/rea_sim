args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

subdir = paste0("sims/", arg.dir, "/", arg.num, "/")

name = "aim2_v2"

## Function to look at data and return posterior probability of being NB, given it is either NB or log-normal
p.nb = function () {
    0.5
}

## Load edgeR and voom-limma with TMM
ed = read.csv(paste0(subdir, "edgeR_res.csv"), row.names=1)
vd = read.csv(paste0(subdir, "voom_TMM_res.csv"), row.names=1)

genes = intersect(rownames(ed), rownames(vd))
ed = ed[genes,]
vd = vd[genes,]

## Average test statistics
pnb = p.nb()
new.log2fc = pnb*ed$logFC + (1-pnb)*vd$log2FC
new.p = pnb*ed$PValue + (1-pnb)*vd$p.value

nd = data.frame(log2FC=new.log2fc, p.value=new.p)
rownames(nd) = genes
nd = nd[order(nd$p.value),]

write.csv(nd, paste0(subdir, name, "_res.csv"))
