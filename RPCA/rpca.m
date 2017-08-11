%% Robust PCA 

addpath(fullfile('~','inexact_alm_rpca'));
addpath(fullfile('~','inexact_alm_rpca', 'PROPACK'));
addpath(fullfile('~','emily_functions','projection'))
addpath(fullfile('~','emily_functions'))

% input and gene cutoff (if necessary) 
inputF = importdata('GSM1599497_ES_d2_LIFminus.csv');
data_full = inputF.data;

%% normalization - choose the desired data pre-processing scheme
%library size  
libsize  = sum(transpose(data_full),2);
data_norm = bsxfun(@rdivide, transpose(data_full), libsize) * median(libsize);
data_norm = data_norm';

% CLR - run in R 
% data_norm = importdata('gsm1599497_clr.txt');

% z-score 
% data_norm = zscore(data_full')';
% 
% log2(FC)
% dataMean = mean(data_full+1,2);
% data_norm = log2((data_full+1)./repmat(dataMean,1,numCells));
    
%% rpca 
[numGenes, numCells] = size(data_norm);
lambda = 1/sqrt(max(numGenes, numCells));
tol = 1e-7;
maxIter = 1000;
[A_hat E_hat iter] = inexact_alm_rpca(data_norm, lambda, tol, maxIter);

save('gsm97_rpca.mat', 'A_hat'); 

%% Generate heatmap

% z-score 
zA_hat_row = zscore(A_hat')';
zA_hat_col = zscore(A_hat');
pcs = 4; % number of desired clusters 
A_hat_plot = zA_hat_row;

pcData = zeros(numGenes, numCells);

pcData_row = zA_hat_row;
pcData_col = zA_hat_col;

% k-means

kids_rows = kmeans(pcData_row,pcs);
kids_cols = kmeans(pcData_col, pcs);

% cluster within clusters - row direction 
inds_ordered_rows = zeros(numGenes,1);
ex_ids_rows = unique(kids_rows);
clust_num_rows = length(ex_ids_rows);
start_spots_rows = zeros(clust_num_rows,1);
clust_sizes_rows = zeros(clust_num_rows,1);
start_rows = 1;
figure(3), clf
for clust = 1:clust_num_rows
    clust_ind_rows = find(kids_rows == ex_ids_rows(clust));
    if length(clust_ind_rows) > 1
        pdis = pdist(pcData_row(clust_ind_rows,:));
        link = linkage(pdis,'average');
        [h t horder] = dendrogram(link,0);
    else
        horder = 1;
    end
    clust_size_rows = length(clust_ind_rows);
    clust_sizes_rows(clust) = clust_size_rows;
    % order indices
    inds_ordered_rows(start_rows:start_rows+clust_size_rows-1) = clust_ind_rows(horder);
    start_rows = clust_size_rows+start_rows;
    start_spots_rows(clust) = start_rows;
end

% cluster within clusters - column direction 
inds_ordered_cols = zeros(numCells,1);
ex_ids_cols = unique(kids_cols);
clust_num_cols = length(ex_ids_cols);
start_spots_cols = zeros(clust_num_cols,1);
clust_sizes_cols = zeros(clust_num_cols,1);
start_cols = 1;
figure(4), clf
for clust = 1:clust_num_cols
    clust_ind_cols = find(kids_cols == ex_ids_cols(clust));
    if length(clust_ind_cols) > 1
        pdis = pdist(pcData_col(clust_ind_cols,:));
        link = linkage(pdis,'average');
        [h t horder] = dendrogram(link,0);
    else
        horder = 1;
    end
    clust_size_cols = length(clust_ind_cols);
    clust_sizes_cols(clust) = clust_size_cols;
    % order indices
    inds_ordered_cols(start_cols:start_cols+clust_size_cols-1) = clust_ind_cols(horder);
    start_cols = clust_size_cols+start_cols;
    start_spots_cols(clust) = start_cols;
end

% plot clustering results 
figure (5), clf 
imagesc(A_hat_plot(inds_ordered_rows,inds_ordered_cols))
colormap redblue
colorbar
ax = axis();
hold on
for clust = 1:clust_num_rows
    plot([start_spots_rows(clust) start_spots_rows(clust)]-.5,[ax(3) ax(4)],'k',...
        'LineWidth',2)
end
hold on
for clust = 1:clust_num_cols
    plot([start_spots_cols(clust) start_spots_cols(clust)]-.5,[ax(3) ax(4)],'k',...
        'LineWidth',2)
end
title(['Expression of K-means ' num2str(pcs) ' clusters'],'Fontsize',14)

