args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

th = c(0.01, 0.05, 0.20)

checker = function (name) {
    de = read.csv(paste0("sims/", arg.dir, "/", arg.num, "/", name,"_res.csv"), row.names=1)
    nam = names(de)
    
    if ("logCPM" %in% nam) {
        ## edgeR
        all = de$PValue
        expressed = de$logCPM > log2(3.5+2)
    } else if ("lfcSE" %in% nam) {
        ## DESeq2
        all = de$pvalue
        expressed = de$baseMean > 3.5
    } else if ("mean" %in% nam) {
        ## t-test
        all = de$p.value
        expressed = de$mean > 3.5
    } else {
        ## Other = voom
        all = de$p.value
        expressed = de$baseMean > log2(3.5+0.505)
    }
    exp = all
    exp[is.na(expressed) | ! expressed] = NA
    
    n.exp = sum(expressed, na.rm=TRUE)
    n.all = nrow(de)

    ## Calc adjusted pvals
    exp.fdr = p.adjust(exp, method="BH")
    all.fdr = p.adjust(all, method="BH")
    exp.bon = p.adjust(exp, method="bonferroni")
    all.bon = p.adjust(all, method="bonferroni")

    ## Load metadata showing which genes had effects induced
    rows = read.table(paste0("sims/", arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
    dex.genes = rownames(rows)[rows$log2FC != 0]
    dex.exp = intersect(dex.genes, rownames(de)[expressed])

    dex.ind = rownames(de) %in% dex.genes
    
    ## Start generating output
    n.dex = length(dex.genes)
    n.dex.exp = length(dex.exp)
    ## out = list(n.exp=n.exp, dex.exp=n.dex.exp)
    out = c(n.all.exp=n.exp, n.dex.exp=n.dex.exp)

    ## Generate output by threshold for 1 analysis
    for (t in th) {
        ## Count sig ----------------------------------------
        ## All genes
        s.un = sum(all < t, na.rm=TRUE)
        s.fdr = sum(all.fdr < t, na.rm=TRUE)
        s.bon = sum(all.bon < t, na.rm=TRUE)

        ## Exp genes
        se.un = sum(exp < t, na.rm=TRUE)
        se.fdr = sum(exp.fdr < t, na.rm=TRUE)
        se.bon = sum(exp.bon < t, na.rm=TRUE)

        ## Count sig dex -------------------------------------
        ## All genes
        sd.un = sum(all[dex.ind] < t, na.rm=TRUE)
        sd.fdr = sum(all.fdr[dex.ind] < t, na.rm=TRUE)
        sd.bon = sum(all.bon[dex.ind] < t, na.rm=TRUE)

        ## Exp genes
        sde.un = sum(exp[dex.ind] < t, na.rm=TRUE)
        sde.fdr = sum(exp.fdr[dex.ind] < t, na.rm=TRUE)
        sde.bon = sum(exp.bon[dex.ind] < t, na.rm=TRUE)

        ## Add to output
        nam = names(out)
        out = c(out,
                s.un, s.fdr, s.bon,
                se.un, se.fdr, se.bon,
                sd.un, sd.fdr, sd.bon, 
                sde.un, sde.fdr, sde.bon)
        names(out) = c(nam, paste0(c("s.un", "s.fdr", "s.bon",
                                     "se.un", "se.fdr", "se.bon",
                                     "sd.un", "sd.fdr", "sd.bon",
                                     "sde.un", "sde.fdr", "sde.bon"), sprintf("%0.2f", t)))
    }

    ## Output
    out
}

modes = c("deseq2", "edgeR", "voom", "ttest", "deseq2_notrim")

out = sapply(modes, checker)
write.table(out, file=paste0("sims/", arg.dir, "/", arg.num, "/check2.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
