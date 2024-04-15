# Generalized Matrix Decomposition Framework (GMDF)


GMDF is a method for unsupervised meta-analysis of diverse datasets, including single-cell genomics data.

**Guidelines:**

See GMDF_stepByStepExample.R script to set up the environment, generate pan-cancer shared programs (using N1 = 25), and run a toy GMDF example. 

To reproduce the CD8 pan-cancer programs: 
source(“GMDF_wrapper.R")
rslts<-GMDF_combine_pancancer()

To use GMDF for other datasets:
rslts<-GMDF_wrapper(E, a, k, k1,N1 ,outputdir)

Annotations:
•	g = number of genes
•	k = number of shared programs
•	k1 = number of context-specific programs
•	n = number of datasets
•	m = number of contexts (e.g, cancer types)

Input:
•	E: list of n gene expression matrices, one per dataset/condition. 
•	a: n x m matrix denoting the value of m annotations (i.e., contexts) in n datasets provided in E (see examples below).**
•	k: Initial estimate of the number of shared programs.
•	K1: Initial estimate of number of context-specific programs per context.
•	N1: Number of times to run GMDF. If N1 > 1 then multiple solutions will be combined.
•	outputdir: The output directory to save the single GMDF run results. Required if N1 > 1.
 
Output description for N1 > 1, the results of N1 combined GMDF runs:
•	sig – the final shared signatures (top 100 genes in the programs represented in the Wf matrix).
•	Wf – the final GMDF shared programs based on all the GMDF runs, such that Wf[i,j] denotes the weight of gene i in program j.
•	W.multi.run – weight of genes (rows) in each of the shared programs identified in the single GMDF runs (columns)
•	W.clusters – the clustering annotation of each one of the shared GMDF programs from the single runs

Output description for N1 = 1, the results of a single GMDF run:
The decomposition of the data
Et[[i]] ~ (Hw[[i]]) x W + ∑j(a[i,j] x  H[[j]][[i]] x A[[j]])
Where Et is transpose(E[[i]][g,]). For additional information about the output and its interpretation see Figure 1A and equation (1) as provided in the STAR METHODS describing GMDF.
  
Output when N1 = 1:
•	Hw - shared programs usage per dataset: list of n matrix, one per dataset, each of size g x k
•	W - g x k matrix representing the shared programs.
•	A - context-specific programs, k1 x g matrix per each of the m contexts
•	H[[j]] - context-specific program usage per context j, with a (g x k1) matrix per each of the n datasets
•	obj - the final value of the objective function, minimizing the reconstruction error
•	objAA - the values of the objective function in each iteration; NA in case the run terminated before the 100th iteration.
•	param - the parameters that were provided as input and used
•	input.data - the input object

**Dependencies:** plyr, rliger. 
