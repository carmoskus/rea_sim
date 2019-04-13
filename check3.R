args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.num = args[2]

if (is.na(arg.dir) || is.na(arg.num) || nchar(arg.dir) == 0 || nchar(arg.num) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

th = c(2^-(200:1), 1-2^-(2:50))

checker = function (name) {
    de = read.csv(paste0("sims/", arg.dir, "/", arg.num, "/", name,"_res.csv"), row.names=1)
    nam = names(de)
    
    if ("logCPM" %in% nam) {
        ## edgeR
        all = de$PValue
    } else if ("lfcSE" %in% nam) {
        ## DESeq2
        all = de$pvalue
    } else if ("mean" %in% nam) {
        ## t-test
        all = de$p.value
    } else {
        ## Other = voom
        all = de$p.value
    }
    
    ## Load metadata showing which genes had effects induced
    rows = read.table(paste0("sims/", arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
    dex.genes = rownames(rows)[rows$log2FC != 0]
    nonzero.genes = rownames(rows)[!is.na(rows$mean) & rows$mean != 0]
    exp.genes = rownames(rows)[!is.na(rows$mean) & rows$mean >= 3.5]
    
    dex.ind = rownames(de)[rownames(de) %in% nonzero.genes] %in% dex.genes

    ## Remove mean=zero genes from all
    all = all[rownames(de) %in% nonzero.genes]
    
    ## Setup list of expressed pvalues
    expressed = rownames(de)[rownames(de) %in% nonzero.genes] %in% exp.genes
    n.exp = sum(expressed, na.rm=TRUE)
    n.all = nrow(de)
    exp = all
    exp[is.na(expressed) | ! expressed] = NA
    dex.exp = intersect(dex.genes, exp.genes)

    ## Calc adjusted pvals
    exp.fdr = p.adjust(exp, method="BH")
    all.fdr = p.adjust(all, method="BH")
    exp.bon = p.adjust(exp, method="bonferroni")
    all.bon = p.adjust(all, method="bonferroni")

    ## Start generating output
    n.dex = length(dex.genes)
    n.dex.exp = length(dex.exp)
    n.nonzero = length(nonzero.genes)
    n.dex.nonzero = length(intersect(dex.genes, nonzero.genes))
    
    out = c(n.all.exp=n.exp, n.nonzero=n.nonzero, n.dex.exp=n.dex.exp, n.dex=n.dex, n.dex.nonzero=n.dex.nonzero)

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
                                     "sde.un", "sde.fdr", "sde.bon"), sprintf("%a", t)))
    }

    ## Output
    out
}

## modes = c("deseq2", "edgeR", "voom", "ttest", "deseq2_notrim", "voom_quantile", "voom_TMM", "ttest_none", "deseq2_ns", "ttest_log")
## modes = c("deseq2", "edgeR", "voom", "ttest", "deseq2_notrim", "voom_quantile", "voom_TMM", "ttest_log")
## modes = c("deseq2", "deseq2_ns", "deseq2_trim_ms",
##           "deseq2_notrim", "deseq2_notrim_ms", "deseq2_notrim_ns",
##           "edgeR",
##           "voom", "voom_ns", "voom_quantile", "voom_TMM",
##           "ttest", "ttest_none",
##           "ttest_log", "ttest_log_ns", "ttest_log_TMM",
##           "ttest_vscale", "ttest_vnone", "ttest_vquantile",
##           "aldex_all_1", "aldex_all_2", "aldex_iqlr_1", "aldex_iqlr_2")
## modes = c("deseq2_notrim", "deseq2_notrim_ms",
##           "edgeR","edgeR_ms",
##           "voom", "voom_TMM",
##           "ttest_log", "ttest_log_TMM",
##           "aldex_all_1", "aldex_all_2", "aldex_iqlr_1", "aldex_iqlr_2")
modes = c("deseq2_notrim",
          "edgeR", "edgeR_notrend",
          "voom_TMM",
          "ttest_log_TMM",
          "aldex_iqlr_1", "aldex_iqlr_2")

out = sapply(modes, checker)
write.table(out, file=paste0("sims/", arg.dir, "/", arg.num, "/check3.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
