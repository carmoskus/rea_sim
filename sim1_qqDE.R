
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
        adj.pvals3 = p.adjust(exp.pvals, method="bonferroni")
        adj.pvals4 = p.adjust(pvals, method="bonferroni")
        
        th = c(0.01, 0.05, 0.20)
        # Output counts of failed tests and lambdas
        c(exp.lambda, n / n2 * 100, all.lambda,
          n, n2,
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
    names(r) = c("Exp at 0.5", "Exp at 0.75", "Exp at 0.9", "Exp at 0.95", "Exp at 0.99",
             "% Exp",
             "All at 0.5", "All at 0.75", "All at 0.9", "All at 0.95", "All at 0.99",
             "Num Expressed", "Num Total",
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
r2 = analyze(range, "edgeR")
r3 = analyze(range, "voom")
r4 = analyze(range, "ttest")
r5 = analyze(range, "deseq2_notrim")

Sys.sleep(3)

colMeans(r1)
r = r1
for (i in 1:11)
    print(c(min(r[,i]), median(r[,i]), max(r[,i])))

colMeans(r2)
r = r2
for (i in 1:11)
    print(c(min(r[,i]), median(r[,i]), max(r[,i])))

colMeans(r3)
r = r3
for (i in 1:11)
    print(c(min(r[,i]), median(r[,i]), max(r[,i])))

colMeans(r4)
r = r4
for (i in 1:11)
    print(c(min(r[,i]), median(r[,i]), max(r[,i])))

colMeans(r5)
r = r5
for (i in 1:11)
    print(c(min(r[,i]), median(r[,i]), max(r[,i])))

data.frame(colnames(r1))
colMeans(r1)
colMeans(r1[r1[,14] > 0,])

colMeans(r1[r1[,15] > 0,])

# Analyses with more DEX genes have higher inflation
colMeans(r1)
colMeans(r1[r1[,16] > 0,])
colMeans(r1[r1[,16] > 1,])

# Analyses with higher inflation have more DEX genes
colMeans(r1)
colMeans(r1[r1[,5] > 1,])
colMeans(r1[r1[,5] < 1,])
