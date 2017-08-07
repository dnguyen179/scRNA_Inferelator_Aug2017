# scRNA_Inferelator_Aug2017
This repository contains source code for project "Transcriptional regulatory network inferrence from single-cell RNA measurements of embryonic stem cells" 

It includes:
1. Plot genes vs percent of zeros 
2. Principle Component Analysis 
3. MAGIC (Markov Affinity-based Graph Imputation of Cells 
4. BISCUIT (Bayesian Inference for Single-cell Clustering and Imputing)
5. RPCA (Robust Principal Component Analysis)
6. Random split gold-standard (G.S) matrix into G.S and prior matrix 
7. Post-Inferelator

I. Preliminary data visualization and principle component analysis
- Plot number of genes according to percent of zeros to determine cutoff threshold for subsequent analysis. This ensures that selected genes are expressed in the dataset (gene_vs_perZero.m; plot_gene_zero.m)
- Use principle component analysis (scores and loadings plots) to show how much of the trends can be explained by principle components. Scores plots are visualized by percent of zeros  (pca_anal_example.m)
+ input: data matrix with gene names and conditions (if applicable) 

II. Data normalization and imputation methods 
1. MAGIC (code: magic.m)
a. Running the package
- Set parameters: 
+ npca: number of principle components 
+ ka: number of nearest neighbors 
+ k : minimum number of neighbors for each cell 
+ t: diffusion time, usually between 6 – 12, smaller ka/k requires bigger t 
+ library size normalization: TRUE or FALSE (normally TRUE)
+ log transformation: TRUE or FALSE (normally FALSE)

- Input data matrix has to be cells x genes 

- Function run_magic(data, t, npca, ka, k, lib_size_norm, log_transform) returns imputed data matrix 

b. Scaling data 
- Calculate the median and mean/median absolute deviation (MAD) for nonzero values in original data (normalized by median of library size) and imputed data from MAGIC
- Scaling factors are the difference of medians and the ratio of MADs
- Scaled data = (MAGIC-imputed data + difference of medians) * ration of MADs 

c. Data visualization 
- Distribution of log gene expression of original data and MAGIC-imputed data (+ pseudocount) via histogram 
- Density plot of original data and imputed data to show the correlation
- Heatmap of cluster-specific gene expression 
+ Need to change the number of desired clusters, default = 4 
+ z-score across genes and samples 
+ using k-means clustering for both row and column direction 

2. RPCA (rpca.m)
a. Pre-processing data 
- by library size (median) 
- by centered log-ration transformation (CLR) (R package)
+ install package “scone” in R
+ use function CLR_FN to return a CLR-normalized data matrix 

b. Running the package 
- Download from: 
http://perception.csl.illinois.edu/matrix-rank/sample_code.html
- Use inexact augmented lagrange multiplier (ALM) method
- Calculate lamda = 1/sqrt(max(number of genes, number of cells))
- Function inexact_alm_rpca(data matrix, lamda, tolerance for stopping criterion, maximum number of iterations) returns estimates of matrix A (A_hat) and E (E_hat)

c. Data visualization
- Heatmap of cluster-specific gene expression 
+ z-score across genes and samples 
+ using k-means clustering for both row and column direction 

3. BISCUIT
- Download BISCUIT package on: https://github.com/sandhya212/BISCUIT_SingleCell_IMM_ICML_2016

- In R:
+ install packages: install.packages(c("MCMCpack","mvtnorm","ellipse","coda","Matrix","Rtsne","gtools","foreach","doParallel","doSNOW","snow","lattice","MASS","bayesm","robustbase","chron","mnormt","schoolmath","devtools","RColorBrewer"))
+ set working_path = getwd()

- Edit the following files:
+ start_file.R: input file name, number of cells, number of genes, number of genes per batch, number of iterations, number of cores, labels of cells (TRUE or FALSE), number of cells per batch
+ BISCUIT_process_data.R: 
•	Change HEADER = FALSE (line 54 or 56 depending on input file type)
•	Choose 1 of the 3 methods to get meaningful genes (usually choose Idea 1 based on standard deviation)
- Output:
+ export alpha_inferred_final, beta_inferred_final, z_inferred_final 
+ export selected gene expression data matrix (when perform one of the three methods): X_all (already log-transformed)
+ export imputed data matrix: Y_rt_final 

- BISCUIT calculations (BISCUIT_calculations.m)
+ Impute data matrix based on linear transformation y = Ax+b 
+ Generate heatmap based on clusters inferred from BISCUIT (z_inferred_final) and density plot to show correlation


III. Randomly split 50% edges of G.S into 1 G.S and 1 prior 
rand_GS_50.m
- For each TF, find the number of edges (non-zero values) 
- Split 50% of the edges to G.S matrix and the other 50% to prior matrix
- Output: 2 files – GS and prior 

IV. Post-Inferelator
- Expand prediction groups into transcription factors that belong to each group with respective edges in the network 
