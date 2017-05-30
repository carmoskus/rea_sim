args = commandArgs(trailingOnly=TRUE)
arg1 = as.integer(args[1])
arg2 = as.integer(args[2])

checker = function (name) {
    function (i) {
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
        
        adj.pvals1 = p.adjust(exp.pvals, method="BH")
        adj.pvals2 = p.adjust(pvals, method="BH")
        adj.pvals3 = p.adjust(exp.pvals, method="bonferroni")
        adj.pvals4 = p.adjust(pvals, method="bonferroni")
        
        n = sum(expressed, na.rm=TRUE)
        n2 = nrow(de)

        th = c(0.01, 0.05, 0.20)
        # Output counts of failed tests
        c(sum(adj.pvals1 < th[1], na.rm=TRUE),
          sum(adj.pvals1 < th[2], na.rm=TRUE),
          sum(adj.pvals1 < th[3], na.rm=TRUE),
          sum(adj.pvals2 < th[1], na.rm=TRUE),
          sum(adj.pvals2 < th[2], na.rm=TRUE),
          sum(adj.pvals2 < th[3], na.rm=TRUE),
          n, n2,
          sum(adj.pvals3 < th[1], na.rm=TRUE),
          sum(adj.pvals3 < th[2], na.rm=TRUE),
          sum(adj.pvals3 < th[3], na.rm=TRUE),
          sum(adj.pvals4 < th[1], na.rm=TRUE),
          sum(adj.pvals4 < th[2], na.rm=TRUE),
          sum(adj.pvals4 < th[3], na.rm=TRUE))
    }
}

        # Output proportions of failed tests
#        c(sum(adj.pvals1 < th[1], na.rm=TRUE) / n,
#          sum(adj.pvals1 < th[2], na.rm=TRUE) / n,
#          sum(adj.pvals1 < th[3], na.rm=TRUE) / n,
#          sum(adj.pvals2 < th[1], na.rm=TRUE) / n2,
#          sum(adj.pvals2 < th[2], na.rm=TRUE) / n2,
#          sum(adj.pvals2 < th[3], na.rm=TRUE) / n2,
#          n, n2,
#          sum(adj.pvals3 < th[1], na.rm=TRUE) / n,
#          sum(adj.pvals3 < th[2], na.rm=TRUE) / n,
#          sum(adj.pvals3 < th[3], na.rm=TRUE) / n,
#          sum(adj.pvals4 < th[1], na.rm=TRUE) / n2,
#          sum(adj.pvals4 < th[2], na.rm=TRUE) / n,
#          sum(adj.pvals4 < th[3], na.rm=TRUE) / n)

analyze = function (range, mode) {
    r = sapply(range, checker(mode))
    r = as.data.frame(t(r), stringsAsFactors=FALSE)
    names(r) = c("Exp.FDR < 0.01", "Exp.FDR < 0.05", "Exp.FDR < 0.20",
             "All.FDR < 0.01", "All.FDR < 0.05", "All.FDR < 0.20",
             "Num Expressed", "Num Total",
             "Exp.Bon < 0.01", "Exp.Bon < 0.05", "Exp.Bon < 0.20",
             "All.Bon < 0.01", "All.Bon < 0.05", "All.Bon < 0.20")
    r
}

range = 1:1000

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

# Use t-tests to see if simulated results match expectations
r = r1
t.test(r[[9]], mu=0.01)
t.test(r[[10]], mu=0.05)
t.test(r[[11]], mu=0.20)
t.test(r[[12]], mu=0.01)
t.test(r[[13]], mu=0.05)
t.test(r[[14]], mu=0.20)

# Second version of tests
r = r5
s = matrix(as.integer(r > 0), nrow=nrow(r), ncol=ncol(r))

prop.test(sum(s[,9]), nrow(s), p=0.01)
prop.test(sum(s[,10]), nrow(s), p=0.05)
prop.test(sum(s[,11]), nrow(s), p=0.20)
prop.test(sum(s[,12]), nrow(s), p=0.01)
prop.test(sum(s[,13]), nrow(s), p=0.05)
prop.test(sum(s[,14]), nrow(s), p=0.20)

