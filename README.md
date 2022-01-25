Estimating The Number Differentially Expressed Genes In Two Dependent Experiments

•	Analyzing microarray and RNA-seq data to determine the genes from two dependent experiments that are different in their expressions using R and Excel

•	The pairs of p-values for the two experiments were separately extracted using the student's t-test and moderated t-test

•	Used Histogram-based method to determine the cutoff points for p-values, for comparing the compatible and incompatible expressions and used Intersection method to find differentially expressed genes in both experiments while controlling the false discovery rate (FDR) at different levels using Benjamini-Hochberg for the adjusted p-value

•	Proposed a method using Truncated Bivariate Normal Distribution to produce an estimate that does not depend on FDR to estimate the number of genes that are differentially expressed in the two dependent experiments

•	Used two simulation studies (one involving independent, normally distributed data and one involving microarray data) to compare the performances of our proposed method, and another method proposed in literature to test for consistency of replicate experiments

