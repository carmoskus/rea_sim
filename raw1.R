args = commandArgs(trailingOnly=TRUE)
arg.dir = args[1]

if (is.na(arg.dir) || nchar(arg.dir) == 0) {
    write("Usage: prog.R subdir", stderr())
    quit(save="no", status=1)
}

## Look at raw count data

name = "raw1"
subdir = "sims/"

checker = function (arg.num) {
    ## Load metadata showing which genes had effects induced
    rows = read.table(paste0(subdir, arg.dir, "/", arg.num, "/rows.txt"), row.names=1)
    exp.genes = rownames(rows)[!is.na(rows$mean) & rows$mean >= 3.5]

    ## Load raw counts
    raw = read.table(paste0(subdir, arg.dir, "/", arg.num, "/counts.txt"), row.names=1)
    raw.exp = raw[exp.genes,]

    ## Look at various characteristics
    n.z = sum(rowSums(raw) == 0)
    n.ze = sum(rowSums(raw.exp) == 0)
    n = nrow(raw)
    n.e = nrow(raw.exp)

    c(n, n.e, n.z, n.ze)
}

res = t(sapply(1:1000, checker))
colnames(res) = c("N", "N.E", "N.Z", "N.ZE")

write.table(res, file=paste0(subdir, arg.dir, "/", name, ".txt"))
