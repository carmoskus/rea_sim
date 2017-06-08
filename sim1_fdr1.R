
df = 1

checker = function (name) {
    function (i) {
        subdir = "b"
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
    names(r5) = c("Exp at 0.5", "Exp at 0.75", "Exp at 0.9", "Exp at 0.95", "Exp at 0.99",
             "% Exp",
             "All at 0.5", "All at 0.75", "All at 0.9", "All at 0.95", "All at 0.99",
             
             "Exp.FDR 0.01", "Exp.FDR at 0.05", "Exp.FDR at 0.20",
             "All.FDR 0.01","All.FDR 0.05", "All.FDR 0.20",
             
             "Num Expressed", "Num Total")
    r
}

range = 1:100
range = 1001:2000
range = 2001:3000
range = 3001:4000
range = 4001:5000

range = 1:5000

r1 = analyze(range, "deseq2")
r2 = analyze(range, "edgeR")
r3 = analyze(range, "voom")
r4 = analyze(range, "ttest")
r5 = analyze(range, "deseq2_notrim")

# Do t-tests
i = 1
for (r in list(r1, r2, r3, r4, r5)) {
    print(i)
    print(t.test(r[["Exp.FDR 0.01"]], mu=0.01))
    print(t.test(r[["Exp.FDR at 0.05"]], mu=0.05))
    print(t.test(r[["Exp.FDR at 0.20"]], mu=0.20))
    print(t.test(r[["All.FDR 0.01"]], mu=0.01))
    print(t.test(r[["All.FDR 0.05"]], mu=0.05))
    print(t.test(r[["All.FDR 0.20"]], mu=0.20))
    i = i + 1
}
