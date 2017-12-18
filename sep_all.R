args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]

if (is.na(arg.dir) || nchar(arg.dir) == 0) {
    write("Usage: prog.R subdir num", stderr())
    quit(save="no", status=1)
}

modes = c("deseq2", "edgeR", "voom", "ttest", "deseq2_notrim")
##th = c(0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)
th = c(0.05)

sep.group = function (arg.num) {
    anno = read.table(paste0("sims/",arg.dir,"/",arg.num,"/cols.txt"), row.names=1, header=TRUE, sep="\t")

    sep = function (p) {
        v = mean(p[,1])
        t1 = table(anno$group[p[,1] < v])
        t2 = table(anno$group[p[,1] > v])

        tpr = (max(t1) + max(t2)) / nrow(anno)
        tpr
    }

    l2c = read.csv(paste0("sims/", arg.dir, "/", arg.num, "/deseq2_notrim_log2counts.csv"), row.names=1)
    pca = as.matrix(read.csv(paste0("sims/", arg.dir, "/", arg.num, "/pca1_values.csv"), row.names=1))

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
        
        ## Calc adjusted pvals
        exp.fdr = p.adjust(exp, method="BH")

        de.genes = rownames(de)[!is.na(exp.fdr) & exp.fdr < th]
        if (length(de.genes) == 0) {
            return(NA)
        }
        ##print(de.genes)
        
        p = prcomp(t(l2c[de.genes,]))
        sep(p$x)
    }

    out = c(all=sep(pca))
    out = c(out, sapply(modes, checker))
    out
}

full = t(sapply(1:1000, sep.group))

m = colMeans(full)

write.table(full, file=paste0("sims/", arg.dir, "/sep_all.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(m, file=paste0("sims/", arg.dir, "/sep_means.txt"), row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
