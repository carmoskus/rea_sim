args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

rows = read.table(paste0("sims/", arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
rows = rows[rows$log2FC != 0,]
dex.genes = rownames(rows)
rows.dir = rows$log2FC

if (length(dex.genes) < 2) {
    write("Error: not enough DEX genes", stderr())
    quit(save="no", status=2)
}

checker = function (name) {
    de = read.csv(paste0("sims/", arg.dir, "/", arg.num, "/", name,"_res.csv"), row.names=1)
    nam = names(de)

    de = de[dex.genes,]
    
    if ("logCPM" %in% nam) {
        ## edgeR
        all = de$PValue
        dir = de$logFC
    } else if ("lfcSE" %in% nam) {
        ## DESeq2
        all = de$pvalue
        dir = de$log2FoldChange
    } else if ("mean" %in% nam) {
        ## t-test
        all = de$p.value
        dir = de$log2FC
    } else {
        ## Other = voom
        all = de$p.value
        dir = de$log2FC
    }

    ## Output
    n = nrow(de)
    hit = sum(sign(rows.dir)==sign(dir))
    rmse = sqrt(mean((rows.dir - dir)^2))
    c(n, hit, n-hit, hit/n*100, cor(rows.dir, dir), cor(rows.dir, dir, method="spearman"), rmse)
}

modes = c("edgeR", "voom_TMM", "ttest_log_TMM", "aim2_v6")

out = t(sapply(modes, checker))
colnames(out) = c("N", "Hit", "Miss", "Hit %", "Pearson R", "Spearman R", "RMSE")

write.table(out, file=paste0("sims/", arg.dir, "/", arg.num, "/cor_FC1.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

