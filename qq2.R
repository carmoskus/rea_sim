args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]
arg.ana = args[2]

if (is.na(arg.dir) || is.na(arg.ana) || nchar(arg.dir) == 0 || nchar(arg.ana) == 0) {
    write("Usage: prog.R subdir analysis", stderr())
    quit(save="no", status=1)
}

name = paste0("qq2_", arg.ana)
subdir = "sims/"

checker = function (arg.num) {
    ## Load metadata showing which genes had effects induced
    rows = read.table(paste0(subdir, arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
    exp.genes = rownames(rows)[!is.na(rows$mean) & rows$mean >= 3.5]

    ## Load DEX data and pull out pvalues
    de = read.csv(paste0(subdir,arg.dir,"/",arg.num,"/",arg.ana,"_res.csv"), row.names=1)
    nam = names(de)
    
    if ("logCPM" %in% nam) {
        pvals = de[exp.genes, "PValue"]
    } else if ("lfcSE" %in% nam) {
        pvals = de[exp.genes, "pvalue"]
    } else if ("mean" %in% nam) {
        pvals = de[exp.genes, "p.value"]
    } else {
        pvals = de[exp.genes, "p.value"]
    }

    sort(pvals)
}

pvals = unlist(sapply(1:1000, checker))

write(pvals, file=paste0(subdir,arg.dir,"/",name,".txt"), ncolumns=1)
