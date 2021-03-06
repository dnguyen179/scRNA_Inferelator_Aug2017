# scRNA_network_inference_Aug2017
This repository contains source code for project "Transcriptional regulatory network inference from single-cell RNA measurements of embryonic stem cells"

With the breakthrough of single-cell RNA-seq (scRNA-seq), transcriptome-wide analyses of individual cells can be undertaken, even when sample material is limiting and a priori knowledge of transcription regulation is insufficient. ScRNA-seq quantifies intra-population heterogeneity, and enables study of cell states and gene expression dynamics at a higher resolution, which are masked in bulk RNA-seq measurements. However, current scRNA-seq techniques have low capture rate, with only 5-15% of total transcriptome of each cell. Thus, scRNA-seq data tend to be sparse data matrices with ambiguous zeros, in which sampling zeros (expressed transcripts but not detected, known as “dropout genes”) and true zeros (transcript not expressed) are not well distinguished. For any downstream analysis, it is important to resolve biological variations, the trends of interest, from technical variations.

In this project, we benchmarked four independent statistical approaches to normalize and impute 'dropout genes' in scRNA-seq data matrices: RPM, MAGIC, BISCUIT and RPCA. In order to evaluate the performance of each method, density plots and heatmaps are generated for visualization. Once the dataset has been processed by each of this method, it will be used as the input gene expression matrix for the Inferelator to generate respective network. To assess the recovery of known interactions, we compute the Area Under the Precision Recall (AUPR) and generate Precision-Recall curves to compare with random. The AUPR has a value of 1 when all G.S interactions rank top of the list and close to 0 for random predictions.

--------------

This github repository includes these following files/functions:
1. Plot genes vs percent of zeros 
2. Principle Component Analysis 
3. MAGIC (Markov Affinity-based Graph Imputation of Cells) 
4. BISCUIT (Bayesian Inference for Single-cell Clustering and Imputing)
5. RPCA (Robust Principal Component Analysis)
6. Random split gold-standard (G.S) matrix into G.S and prior matrix 
7. Post-Inferelator

* Software required: MATLAB (tested in R2016a), R (tested in 3.4.0) and Python (tested in 3.6.1)
* All tested data sets are contained in the Datasets folder (in compressed format)
- Uncompress them before use by changing to the

 ../scRNA_network_inferrence_Aug2017-master/Datasets

directory (folder) and issuing this command:

bzip2 -dk *.bz2 <RET>

* For Matlab use, at the start of each program there is a addpath(genpath()) command, which the user will have to alter, substituting in the full path to the downloaded GitHub diretory. For example, if the downloaded directory is in

   /MiraldiLab/scRNA_network_inferrence_Aug2017-master

then the addpath() line should be

addpath(genpath('/MiraldiLab/DiepNguyenGitHubProject/scRNA_network_inferrence_Aug2017-master'))

- Once that is done, the data files should now be ready for use and accessible to the Matlab programs simply by file name.
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
  + To invoke, type this command: 
[in MATLAB] pca_perZeros.m <RET>

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

a. Running the package 
- Download from: 
http://perception.csl.illinois.edu/matrix-rank/sample_code.html
- Use inexact augmented lagrange multiplier (ALM) method
- Change parameters:
  + Input data file 
  + data pre-processing methods: library size, CLR, zscore or log2
  + tol 
  + maxIter 
- Calculate lamda = 1/sqrt(max(number of genes, number of cells))
- Function inexact_alm_rpca(data matrix, lamda, tolerance for stopping criterion, maximum number of iterations) returns estimates of matrix A (A_hat) and E (E_hat)
- Once the input and parameters have been altered, type this command to run the program:

[in MATLAB] rpca.m <RET>

b. Data visualization (done automatically once rpca.m is invoked)
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
    +	Change HEADER = TRUE or FALSE (line 54 or 56, depending on input file type and format)
    +	Choose 1 of the 3 methods to get meaningful genes (usually choose Idea 1 based on standard deviation)
- Output:
  + export alpha_inferred_final, beta_inferred_final, z_inferred_final 
  + export selected gene expression data matrix (when perform one of the three gene selecting methods): X_all (log-transformed)
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
- To invoke, type this command:
[in MATLAB] rand_GS_50.m <RET>

IV. Post-Inferelator
(get_group.R; get_TF.py)
- This step uses R and Python
- Edit the get_group.R R program as desired to select the name of the input Inferelator file. The body of the program will need to be altered to reflect the number of predictor groups. The program code must match the number of groups in the selected Inferelator file used as input.
  + To invoke, type this command: 
[in R]  Rscript get_group.R <RET>
  
- The output of get_group.R is used as the input for get_TF.py: Edit the get_TF.py Python file aS desired to select the name of the input Inferelator file (as the "network"). The "pred_file" will be set to the output file name of the R program run avbove. The "TF_network" will be set to the desired output file name.
  + To invoke, type this command: 
[in Python]  Python3 get_TF.py





