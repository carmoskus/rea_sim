args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.start = as.integer(args[2])
arg.end = as.integer(args[3])

if (is.na(arg.dir) || nchar(arg.dir) == 0 || is.na(arg.start) || is.na(arg.end)) {
    write("Usage: prog.R subdir num.start num.end", stderr())
    quit(save="no", status=1)
}

th = c(2^-(200:6), 2^-6*2:62, 1-2^-(6:50))

checker = function (arg.dir) function (name) {
    de = read.csv(paste0("sims/", arg.dir, "/", arg.num, "/", name,"_res.csv"), row.names=1)
    nam = names(de)
    
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
    
    ## Load metadata showing which genes had effects induced
    rows = read.table(paste0("sims/", arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
    dex.genes = rownames(rows)[rows$log2FC != 0]
    nonzero.genes = rownames(rows)[!is.na(rows$mean) & rows$mean != 0]
    exp.genes = rownames(rows)[!is.na(rows$mean) & rows$mean >= 3.5]
    
    dex.ind = rownames(de)[rownames(de) %in% nonzero.genes] %in% dex.genes

    ## Remove mean=zero genes from all
    all = all[rownames(de) %in% nonzero.genes]
    dir = dir[rownames(de) %in% nonzero.genes]
    
    ## Setup list of expressed pvalues
    expressed = rownames(de)[rownames(de) %in% nonzero.genes] %in% exp.genes
    n.exp = sum(expressed, na.rm=TRUE)
    n.all = nrow(de)
    exp = all
    exp[is.na(expressed) | ! expressed] = NA
    dex.exp = intersect(dex.genes, exp.genes)
    
    ## Start generating output
    n.dex = length(dex.genes)
    n.dex.exp = length(dex.exp)
    n.nonzero = length(nonzero.genes)
    n.dex.nonzero = length(intersect(dex.genes, nonzero.genes))
    
    out = c(n.all.exp=n.exp, n.nonzero=n.nonzero, n.dex.exp=n.dex.exp, n.dex=n.dex, n.dex.nonzero=n.dex.nonzero)
    
    ## Generate output by threshold for 1 analysis
    de = de[rownames(de) %in% nonzero.genes,]
    for (t in th) {
        ## Count sig ----------------------------------------
        ## All genes
        s.un = sum(all < t, na.rm=TRUE)
        
        ## Exp genes
        se.un = sum(exp < t, na.rm=TRUE)
        
        ## Count sig dex -------------------------------------
        ## All genes
        mask = !is.na(all) & all < t & dex.ind
        ds = rownames(de)[mask]
        dir.real = rows[ds, "log2FC"]
        dir.seen = dir[mask]
        
        sd.un = sum(mask)
        d = (sign(dir.real) == sign(dir.seen))
        sda.un = sum(d, na.rm=TRUE)
        
        ## Exp genes
        b = exp[mask] < t
        sde.un = sum(b, na.rm=TRUE)
        sdae.un = sum(b & d, na.rm=TRUE)
    
        ## Add to output
        nam = names(out)
        out = c(out,
                s.un, se.un, sd.un, sde.un, sda.un, sdae.un)
        names(out) = c(nam, paste0(c("s.un", "se.un", "sd.un", "sde.un", "sda.un", "sdae.un"), sprintf("%a", t)))
    }
    
    ## Output
    out
}

## modes = c("deseq2_notrim", "deseq2_notrim_ms", "deseq2_notrim_ns", "deseq2_notrim_TMM", "deseq2_notrim_UQ", "deseq2_notrim_RLE",
##           "deseq2", "deseq2_trim_ms", "deseq2_ns", 
##           "edgeR", "edgeR_ms", "edgeR_ns", "edgeR_UQ", "edgeR_RLE", "edgeR_notrend",
##           "voom", "voom_ns", "voom_quantile", "voom_TMM", "voom_UQ", "voom_RLE", "voom_manWeights", "voom_onlyE",
##          "ttest", "ttest_none",
##           "ttest_log", "ttest_log_ns", "ttest_log_TMM", "ttest_log_UQ", "ttest_log_RLE", "ttest_vquantile",
##          "ttest_vscale", "ttest_vnone", "ttest_vquantile",
##          "ttest_voomt_TMM", "ttest_voomt_ms", "ttest_voomt_ns",
##          "aim2_v1", "aim2_v2", "aim2_v3", "aim2_v4", "aim2_v5", "aim2_v6") #,
##          "aldex_all_1", "aldex_all_2", "aldex_iqlr_1", "aldex_iqlr_2")
## TODO: update so names include norms as suffix rather than specifying
##s.norms = c("TMM", "RLE", "UQ", "ms", "ns")
## Other norms: DE2 quantile
edgeR.modes = c("edgeR_TMM", "edgeR_RLE", "edgeR_UQ", "edgeR_ms", "edgeR_ns")
voom.modes = c("voom_TMM", "voom_RLE", "voom_UQ", "voom_ms", "voom_ns", "voom_quantile")
d2nt.modes = c("d2notrim_TMM", "d2notrim_RLE", "d2notrim_UQ", "d2notrim_ms", "d2notrim_ns", "d2notrim_DE2")
d2t.modes = c("d2trim_TMM", "d2trim_RLE", "d2trim_UQ", "d2trim_ms", "d2trim_ns", "d2trim_DE2")
ttestL.modes = c("ttest_log_TMM","ttest_log_RLE","ttest_log_UQ","ttest_log_ms","ttest_log_ns")
ttestHL.modes = c("ttestH_log_TMM","ttestH_log_RLE","ttestH_log_UQ","ttestH_log_ms","ttestH_log_ns")
lmFitVT.modes = c("lmFit_TMM","lmFit_RLE","lmFit_UQ","lmFit_ms","lmFit_ns")
#modes = c(edgeR.modes, voom.modes, d2nt.modes, d2t.modes, ttestL.modes, ttestHL.modes, lmFitVT.modes)
modes = c(ttestHL.modes)

lapply(arg.start:arg.end,  function (arg.num) {
    out = sapply(modes, checker(arg.num))
    write.table(out, file=paste0("sims/", arg.dir, "/", arg.num, "/check3.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
})

warnings()
