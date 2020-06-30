args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.start = as.integer(args[2])
arg.end = as.integer(args[3])

if (is.na(arg.dir) || nchar(arg.dir) == 0 || is.na(arg.start) || is.na(arg.end)) {
    write("Usage: prog.R subdir num.start num.end", stderr())
    quit(save="no", status=1)
}

th = c(2^-(200:6), 2^-6*2:62, 1-2^-(6:50))

checker = function (arg.num, name) {
    de = read.csv(paste0("sims/", arg.dir, "/", arg.num, "/", name,"_res.csv"), row.names=1)
    nam = names(de)

    ## Load metadata showing which genes had effects induced
    rows = read.table(paste0("sims/", arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
    dex.genes = rownames(rows)[rows$log2FC != 0]
    nonzero.genes = rownames(rows)[!is.na(rows$mean) & rows$mean != 0]
    exp.genes = rownames(rows)[!is.na(rows$mean) & rows$mean >= 3.5]
    
    dex.ind = rownames(de)[rownames(de) %in% nonzero.genes] %in% dex.genes

    ## Setup list of expressed pvalues
    expressed = rownames(de)[rownames(de) %in% nonzero.genes] %in% exp.genes
    n.exp = sum(expressed, na.rm=TRUE)
    n.all = nrow(de)

    check = function (all, dir) {
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
        out
    }
    
    ## Remove mean=zero genes from all
    all = all[rownames(de) %in% nonzero.genes]
    dir = dir[rownames(de) %in% nonzero.genes]
    exp = all
    exp[is.na(expressed) | ! expressed] = NA
    dex.exp = intersect(dex.genes, exp.genes)
    
    ## Start generating output
    n.dex = length(dex.genes)
    n.dex.exp = length(dex.exp)
    n.nonzero = length(nonzero.genes)
    n.dex.nonzero = length(intersect(dex.genes, nonzero.genes))
        
    out = c(n.all.exp=n.exp, n.nonzero=n.nonzero, n.dex.exp=n.dex.exp, n.dex=n.dex, n.dex.nonzero=n.dex.nonzero)
    
    if ("p13" %in% nam) {
        check(de$p12, z12)
    } else if ("logCPM" %in% nam) {
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
    out
}

modes = c("mix_test1")

o = lapply(arg.start:arg.end,  function (arg.num) {
    out = sapply(modes, function (mode) checker(arg.num, mode))
    write.table(out, file=paste0("sims/", arg.dir, "/", arg.num, "/check_mix.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
})

warnings()
