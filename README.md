rea_sim - Simulation analysis of differential expression testing in RNA-Seq data
========================================================

`rea_sim` helps you to:

- simulate RNA-Seq counts from a variety of statistical models (negative binomial, normal, log-normal, multinomial)
- automate running a variety of differential expression analyses using various programs and combinations of covariates (DESeq2, edgeR, voom-limma, PCA, SVA, PEER, etc.)

Subprograms:

- sim1_create.R: Simulates data from NB distribution. Stores in sims/[a-z]/$n
