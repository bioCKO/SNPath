# SNPath
An R package for pathway analysis of GWAS data (SNPath).

Genome-wide association studies (GWAS) assess the association between individual SNPs and
disease risk, and have successfully identified susceptibility loci for various complex diseases. In
addition, many methods have been proposed to evaluate the association between disease risk and
a set of SNPs that belongs to functional gene sets or pathways. The SNPath package contains four
different algorithms in the literature: grass, gseaSnp, plinkSet and aligator. Users
can use any one of them to identify pathways that are associated with disease risk; meanwhile, this
package provides a nice and convenient platform for comparison of different algorithms as well.


This document provides a tutorial for using the SNPath package, as well as detailed description
of each algorithm in this package. Specifically, the function "grass" calculate the pvalues for a
prior defined pathways (sets of genes), based on the Gene set Ridge regression in ASsociation
Studies (GRASS) algorithm (Chen et al 2010, AJHG, 86(6):860-871). The algorithm summarizes the genetic structure using singular
value decomposition for each gene as eigenSNPs and uses a novel form of regularized regression
technique, termed group ridge regression, to select representative eigenSNPs for each gene and
assess their joint association with disease risk. The function "gseaSnp" implements the phenotype-
permution algorithm proposed in Wang et al. (2007) AJHG 81(6):1278-1283. This algorithm modifies the gene-set
enrichment analysis approach for expression studies and is considered the first approach for gene-
set enrichment analysis in association studies. The function "plinkSet" implements the set-based
tests in the popular whole genome association analysis toolset, PLINK software. The function
aligator performs a simple and fast ALIGATOR algorithm (Holmans et al. AJHG 85(1):12-24), which is a method for testing
overrepresentation of pathways, in lists of significant SNPs from GWAS. With a simulated example,
we will show how to use each function with various options.
