#args = commandArgs(trailingOnly=TRUE)
#arg1 = args[1]

checker = function (name) {
    function (i) {
        de = read.csv(paste0("a/",i,"/",name,"_res.csv"), row.names=1)
        adj.pvals = p.adjust(de$pvalue[de$baseMean >= 3.5], method="BH")
        
        n = sum(de$baseMean >= 3.5)
        n2 = sum(!is.na(de$pvalue))
        
        c(sum(de$padj < 0.01, na.rm=TRUE),
          sum(de$padj < 0.05, na.rm=TRUE),
          sum(de$padj < 0.2, na.rm=TRUE),
          sum(adj.pvals < 0.01),
          sum(adj.pvals < 0.05),
          sum(adj.pvals < 0.20),
            n, n2,
          sum(de$pvalue < 0.01 / n, na.rm=TRUE),
          sum(de$pvalue < 0.05 / n, na.rm=TRUE),
          sum(de$pvalue < 0.20 / n, na.rm=TRUE),
          sum(de$pvalue < 0.01 / n2, na.rm=TRUE),
          sum(de$pvalue < 0.05 / n2, na.rm=TRUE),
          sum(de$pvalue < 0.20 / n2, na.rm=TRUE))
    }
}

range = 1:200
r1 = sapply(range, checker("deseq2"))
rowMeans(r1)
rowSums(r1)
r2 = sapply(range, checker("deseq2_notrim"))
rowMeans(r2)
rowSums(r2)

