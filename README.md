# scRNA_network_inference_Aug2017
This repository contains source code for project "Transcriptional regulatory network inference from single-cell RNA measurements of embryonic stem cells" 

It includes:
1. Plot genes vs percent of zeros 
2. Principle Component Analysis 
3. MAGIC (Markov Affinity-based Graph Imputation of Cells) 
4. BISCUIT (Bayesian Inference for Single-cell Clustering and Imputing)
5. RPCA (Robust Principal Component Analysis)
6. Random split gold-standard (G.S) matrix into G.S and prior matrix 
7. Post-Inferelator

* Software required: MATLAB (tested in R2016a), R (tested in 3.4.0) and Python (tested in 3.6.1)
* All tested data sets are contained in the Datasets folder (in compressed format)

-----------
WORKFLOW

I. Preliminary data visualization and principle component analysis (gene_vs_zeros)
- Plot number of genes according to percent of zeros to determine cutoff threshold for subsequent analysis. This ensures that selected genes are expressed in the dataset (gene_vs_perZero.m; plot_gene_zero.m)
  + Invoke MATLAB
  + Select the data set(s) to visualize (found in ../scRNA_network_inference_Aug2017-master/Datasets)
  + Visualization program: plot_gene_zero.m 
  + Once the desired data file(s) has been hard-coded into the program, to display plot, type this command: 

[in MATLAB] plot_gene_zero.m <RET>

- Principal Component Analysis (in../scRNA_network_inference_Aug2017-master/PCA/pca_perZeros.m): This allows data visualization according to how much of the variability in the data can be accounted by 4 principal components. Data in scores plots are grouped and color-coded by the percent of zeros (ascending by 12.5%) 


II. Data normalization and imputation methods 
1. MAGIC (code: magic.m, in ../scRNA_network_inference_Aug2017-master/MAGIC)

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
- Once the input and parameters have been altered, type this command to run the program:

[in MATLAB] magic.m <RET>

b. Scaling data (done automatically once magic.m is invoked)
- Calculate the median and mean/median absolute deviation (MAD) for nonzero values in original data (normalized by median of library size) and imputed data from MAGIC
- Scaling factors are the difference of medians and the ratio of MADs
- Scaled data = (MAGIC-imputed data + difference of medians) * ration of MADs 

c. Data visualization (done automatically once magic.m is invoked)
- Distribution of log gene expression of original data and MAGIC-imputed data (+ pseudocount) via histogram 
- Density plot of original data and imputed data to show the correlation
- Heatmap of cluster-specific gene expression 
  + Need to change the number of desired clusters, default = 4 
  + z-score across genes and samples 
  + using k-means clustering for both row and column direction 

2. RPCA (rpca.m, in in ../scRNA_network_inference_Aug2017-master/RPCA)

a. Pre-processing data 
- by library size (median) 
- by centered log-ration transformation (CLR) (R package)

b. Running the package 
- Download from: 
http://perception.csl.illinois.edu/matrix-rank/sample_code.html
- Use inexact augmented lagrange multiplier (ALM) method
- Change parameters:
  + Input data file 
  + tol 
  + maxIter 
- Calculate lamda = 1/sqrt(max(number of genes, number of cells))
- Function inexact_alm_rpca(data matrix, lamda, tolerance for stopping criterion, maximum number of iterations) returns estimates of matrix A (A_hat) and E (E_hat)
- Once the input and parameters have been altered, type this command to run the program:

[in MATLAB] rpca.m <RET>

c. Data visualization (done automatically once rpca.m is invoked)
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
  + Input: alpha_inferred_final, beta_inferred_final, z_inferred_final, selected gene expression matrix, mu_final, sigma_final
  + Impute data matrix based on linear transformation y = Ax+b 
  + Generate heatmap based on clusters inferred from BISCUIT (z_inferred_final) and density plot to show correlation
- To run the file, type this command: 

[in MATLAB] BISCUIT_calculations.m <RET>

III. Randomly split 50% edges of G.S into 1 G.S and 1 prior 
(rand_GS_50.m)
- For each TF, find the number of edges (non-zero values) 
- Split 50% of the edges to G.S matrix and the other 50% to prior matrix
- Alter the names of the input file and of the two output files as desired
- Output: 2 files – G.S and prior 
- Once the input and parameters have been altered, type this command to run the program:

[in MATLAB] rand_GS_50.m <RET>

IV. Post-Inferelator
(get_group.R; get_TF.py)
- This step uses R and Python
- Edit the get_group.R R program as desired to select the name of the input Inferelator file. The body of the program will need to be altered to reflect the number of predictor groups. The program code must match the number of groups in the selected Inferelator file used as input.
  
[in R]  Rscript get_group.R <RET>
  
- The output of get_group.R is used as the input for get_TF.py: Edit the get_TF.py Python file aS desired to select the name of the input Inferelator file (as the "network"). The "pred_file" will be set to the output file name of the R program run avbove. The "TF_network" will be set to the desired output file name.
  
[in Python]  Python3 get_TF.py





