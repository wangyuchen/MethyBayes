
<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

DNA methylation is methylation cytosine residues at CpG dinucleotides in a DNA sequence and affects 70-80% of all CpG dinucleotides in mammals. It is the most widely studied epigenetic modification and is known to have profound effects on gene expression. Given the influence of methylation on the gene expression and various type of cancer, there are a lot of studies aiming to identify differentially methylation loci in diseased tissue samples compared to their respective normal samples. DNA methylation change in cancer tissue compared to normal tissue can be increased called hyper-methylation and decreased called hypo-methylation. Since they have very different biological meanings, it is important to differentiate hypo- and hyper-methylation.

Although there are several statistical methods applied to test differential methylated loci, all of them only provide binary outcome of a locus being methylated or not. A common practice to further identify hypo- and hyper-methylated loci is based on the sign of the test statistics, which is an ad-hoc approach and ignores the uncertainty associated with the test statistics. \#\# Implementation

In this MethyBayes package, we implement our new proposed full Bayesian partition model to identify hypo- and hyper-methylated loci simultaneously without further post-hoc analysis. Our proposed Bayesian model has been proved based on the simulation study and real data analysis that it is a both powerful and reliable procedure. The main bayesian routine is implemented in Fortran for efficiency.

Installation
------------

``` r
install.packages("MethyBayes")
```

Example
-------

Here's an example using the included demo datasets:

``` r
library(MethyBayes)
data("case1", "case2", "control1", "control2")
read_counts <- seqdata(list(case1, case2), list(control1, control2))
ident <- MethyBayes(read_counts, 2, 2)
```
