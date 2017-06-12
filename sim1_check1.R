#args = commandArgs(trailingOnly=TRUE)
#arg1 = as.integer(args[1])

checker = function (name) {
    function (i) {
#        de = read.csv(paste0("a/",i,"/",name,"_res.csv"), row.names=1)
        de = read.csv(paste0("a/",i,"/",name,"_res.csv"), row.names=1)
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
        
        adj.pvals3 = p.adjust(exp.pvals, method="bonferroni")
        adj.pvals4 = p.adjust(pvals, method="bonferroni")
        
        n = sum(expressed, na.rm=TRUE)
        n2 = nrow(de)

        th = c(0.01, 0.05, 0.20)
        # Output counts of failed tests
        c(n, n2,
          sum(adj.pvals3 < th[1], na.rm=TRUE),
          sum(adj.pvals3 < th[2], na.rm=TRUE),
          sum(adj.pvals3 < th[3], na.rm=TRUE),
          sum(adj.pvals4 < th[1], na.rm=TRUE),
          sum(adj.pvals4 < th[2], na.rm=TRUE),
          sum(adj.pvals4 < th[3], na.rm=TRUE))
    }
}

analyze = function (range, mode) {
    r = sapply(range, checker(mode))
    r = as.data.frame(t(r), stringsAsFactors=FALSE)
    names(r) = c("Num Expressed", "Num Total",
             "Exp.Bon < 0.01", "Exp.Bon < 0.05", "Exp.Bon < 0.20",
             "All.Bon < 0.01", "All.Bon < 0.05", "All.Bon < 0.20")
    r
}

range = 1:1000
range = 1001:2000
range = 2001:3000
range = 3001:4000
range = 4001:5000

range = 1:5000

r1 = analyze(range, "deseq2")
colMeans(r1)
colMeans(r1 > 0)

r2 = analyze(range, "edgeR")
colMeans(r2)
colMeans(r2 > 0)

r3 = analyze(range, "voom")
colMeans(r3)
colMeans(r3 > 0)

r4 = analyze(range, "ttest")
colMeans(r4)
colMeans(r4 > 0)

r5 = analyze(range, "deseq2_notrim")
colMeans(r5)
colMeans(r5 > 0)

r6 = analyze(range, "d2notrim_sva1")
colMeans(r6)
colMeans(r6 > 0)

r7 = analyze(range, "d2notrim_sva2")
colMeans(r7)
colMeans(r7 > 0)

r8 = analyze(range, "d2notrim_sva3")
colMeans(r8)
colMeans(r8 > 0)

r9 = analyze(range, "d2notrim_sva4")
colMeans(r9)
colMeans(r9 > 0)

r10 = analyze(range, "d2notrim_sva5")
colMeans(r10)
colMeans(r10 > 0)

# Third version of tests
i = 1
for (r in list(r1, r2, r3, r4, r5)) {
    
    s = matrix(as.integer(r > 0), nrow=nrow(r), ncol=ncol(r))

    print(i)
    print(binom.test(sum(s[,3]), nrow(s), p=0.01))
    print(binom.test(sum(s[,4]), nrow(s), p=0.05))
    print(binom.test(sum(s[,5]), nrow(s), p=0.20))
    print(binom.test(sum(s[,6]), nrow(s), p=0.01))
    print(binom.test(sum(s[,7]), nrow(s), p=0.05))
    print(binom.test(sum(s[,8]), nrow(s), p=0.20))

    i = i + 1
}
