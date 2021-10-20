# Generalized Matrix Decomposition Framework (GMDF)


GMDF is a method for unsupervised meta-analysis of diverse datasets, including single-cell genomics data.
To run GMDF use call_GMDF with the following input:

- L (nx1) = list of n gene expression matrixes, one per dataset
- a (nxk) = matrix of convariates, a(i,j) denoting covariate j in dataset i
- n.shared (default 5): pre-defined number of shared programs
- n.spc (default 5): pre-defined number of covariate-specific programs
- var.thresh (default 0.3): Threshold used to identify variable genes (higher threshold -> fewer selected genes). 
- output.file (optional): The file where the results will be saved
