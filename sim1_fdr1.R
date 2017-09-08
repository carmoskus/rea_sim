
df = 1
subdir = "b"

checker = function (name) {
    function (i) {
        de = read.csv(paste0(subdir, "/",i,"/",name,"_res.csv"), row.names=1)
        nam = names(de)
        
        if ("logCPM" %in% nam) {
            # edgeR
            pvals = de$PValue
            expressed = de$logCPM > log2(3.5+2)
        } else if ("lfcSE" %in% nam) {
            # DESeq2
            pvals = de$pvalue
            expressed = de$baseMean > 3.5
        } else if ("mean" %in% nam) {
            # t-test
            pvals = de$p.value
            expressed = de$mean > 3.5
        } else {
            # Other = voom
            pvals = de$p.value
            expressed = de$baseMean > log2(3.5+0.505)
        }
        exp.pvals = pvals
        exp.pvals[is.na(expressed) | ! expressed] = NA
        
        n = sum(expressed, na.rm=TRUE)
        n2 = nrow(de)

        quantiles = c(0.5, 0.25, 0.1, 0.05, 0.01)
        
        # Calculate stats for expressed
        exp.p = quantile(exp.pvals, quantiles, na.rm=TRUE)
        exp.chi = qchisq(1-exp.p, df)
        exp.lambda = exp.chi / qchisq(1-quantiles, df)

        # Calculate stats for all
        all.p = quantile(pvals, quantiles, na.rm=TRUE)
        all.chi = qchisq(1-all.p, df)
        all.lambda = all.chi / qchisq(1-quantiles, df)

        #  Calc adjusted pvals
        de$adj.pvals1 = adj.pvals1 = p.adjust(exp.pvals, method="BH")
        de$adj.pvals2 = adj.pvals2 = p.adjust(pvals, method="BH")

        # Load metadata showing which genes had effects induced
        rows = read.table(paste0(subdir, "/",i,"/rows.txt"), row.names=1)
        dex.genes = rownames(rows)[rows$log2FC != 0]
        dex.exp = intersect(dex.genes, rownames(de)[expressed])

        DEX = 0
        de[dex.genes,"DEX"] = 1
        
        th = c(0.01, 0.05, 0.20)
        d1 = c()
        d2 = c()
        td1 = c()
        td2 = c()
        for (i in 1:length(th)) {
            t = th[i]
            d1[i] = sum(adj.pvals1 < th[i], na.rm=TRUE)
            d2[i] = sum(adj.pvals2 < th[i], na.rm=TRUE)
            td1[i] = sum(de[dex.genes, "adj.pvals1"] < th[i], na.rm=TRUE)
            td2[i] = sum(de[dex.genes, "adj.pvals2"] < th[i], na.rm=TRUE)
        }
        fd1 = d1 - td1
        fd2 = d2 - td2

        fdr1 = fd1 / d1
        fdr2 = fd2 / d2
        
        # Output counts of correct and false tests and lambdas
        c(exp.lambda, n / n2 * 100, all.lambda,
          fdr1, fdr2,
          n, n2)
    }
}

analyze = function (range, mode) {
    r = sapply(range, checker(mode))
    r = as.data.frame(t(r), stringsAsFactors=FALSE)
    names(r) = c("Exp at 0.5", "Exp at 0.75", "Exp at 0.9", "Exp at 0.95", "Exp at 0.99",
             "% Exp",
             "All at 0.5", "All at 0.75", "All at 0.9", "All at 0.95", "All at 0.99",
             
             "Exp.FDR 0.01", "Exp.FDR at 0.05", "Exp.FDR at 0.20",
             "All.FDR 0.01","All.FDR 0.05", "All.FDR 0.20",
             
             "Num Expressed", "Num Total")
    r
}

range = 1:5000

modes = c("deseq2", "edgeR", "voom", "ttest", "deseq2_notrim")

modes = c(modes, paste0("d2notrim_pca1n", 1:5))
modes = c(modes, paste0("d2notrim_sva", 1:5))
modes = c(modes, paste0("edgeR_pca1", 1:5))
modes = c(modes, paste0("edgeR_sva", 1:5))

args = commandArgs(trailingOnly=TRUE)
arg1 = as.integer(args[1])

if (!is.na(arg1) || arg1 < 1) {
    r = analyze(range, modes[arg1])
    # Output results
    write.table(r, file=paste0("summaries/", subdir, "1-5k_FDR_", modes[arg1], ".txt"), sep="\t", quote=FALSE)
}

quit(save="no")

r1 = analyze(range, "deseq2")
r2 = analyze(range, "edgeR")
r3 = analyze(range, "voom")
r4 = analyze(range, "ttest")

r5 = analyze(range, "deseq2_notrim")

r6 = analyze(range, "d2notrim_sva1")
r7 = analyze(range, "d2notrim_sva2")
r8 = analyze(range, "d2notrim_sva3")
r9 = analyze(range, "d2notrim_sva4")
r10 = analyze(range, "d2notrim_sva5")

# Do t-tests
Sys.sleep(2)
i = 0
#for (r in list(r1, r2, r3, r4, r5)) {
for (r in list(r5, r6, r7, r8, r9, r10)) {
    print(i)
    print(t.test(r[["Exp.FDR 0.01"]], mu=0.01))
    print(t.test(r[["Exp.FDR at 0.05"]], mu=0.05))
    print(t.test(r[["Exp.FDR at 0.20"]], mu=0.20))
    print(t.test(r[["All.FDR 0.01"]], mu=0.01))
    print(t.test(r[["All.FDR 0.05"]], mu=0.05))
    print(t.test(r[["All.FDR 0.20"]], mu=0.20))
    i = i + 1
}
