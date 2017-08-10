%% BISCUIT calculations 

addpath(fullfile('~','Desktop', 'BISCUIT_aug7','output', 'plots', 'Inferred_Sigmas'))
addpath(fullfile('~','Desktop', 'BISCUIT_aug7','output', 'plots', 'Inferred_means'))
addpath(fullfile('~','Desktop', 'BISCUIT_aug7'));
addpath(fullfile('~','Documents','emily_functions','projection'))
addpath(fullfile('~','Documents','emily_functions'))


cluster_num = importdata('z_inferred_final.txt');
data_raw = importdata('selected.txt');
alpha_full = importdata('alpha_inferred_final.txt');
beta_full = importdata('beta_inferred_final.txt');
mean_full = importdata('mu_final.txt'); 
data_full = transpose(data_raw);

[numGenes, numCells] = size(data_full);
[var, totK] = size(mean_full);
newData = zeros(numGenes, numCells);

% data_i = importdata('imputed.txt');
% data_imputed = data_i';

for i = 1:numCells
    beta_cell = beta_full(i);
    alpha_cell = alpha_full(i);
    cluster = cluster_num(i);
    FileName = strcat('Sigma_final_', num2str(cluster),'.txt');
    sigma = importdata(FileName);
    sigma_new = (1/beta_cell)*sigma;
    [v_cell, s_cell, null] = svd(sigma_new);
    A_cell = v_cell*(s_cell)^(1/2)*((sigma)^(-1/2));
    b_cell = mean_full(:,cluster) - alpha_cell*A_cell*mean_full(:,cluster);
    x = data_full(:,i);
    newData(:,i) = A_cell*x + b_cell;

end


%% Categorize and generate heatmap

% imputed = importdata('imputed.txt');
% data_imputed = imputed';


ClusterMatrix = [];
for k = 1:totK
    for i = 1: numCells
        if cluster_num(i) == k 
            ClusterMatrix = [ClusterMatrix, data_imputed(:,i)];
        end
    end
end

zClusterMatrix = zscore(ClusterMatrix')';
geneDistExp = pdist(zClusterMatrix);
geneLinkExp = linkage(geneDistExp, 'ward');
[h, t, horderGeneExp] = dendrogram(geneLinkExp, 0);
figure(22)
subplot(1,5,2:5)
colormap redblue
imagesc(zClusterMatrix(horderGeneExp,:))
title('expression')


for i = 1:numGenes
    currExp = zscore(cluster_final(:,i)')';
    rawExp = cluster_final(:,i);

    geneDistExp = pdist(zClusterMatrix);
    geneLinkExp = linkage(geneDistExp, 'ward');
    [h, t, horderGeneExp] = dendrogram(geneLinkExp, 0);
    figure(22)
    subplot(1,5,2:5)
    colormap redblue
    imagesc(zClusterMatrix(horderGeneExp,:))
    title('expression')
end

exp = [];
for i = 1:numGenes
    aveExp = zscore(cluster_final(:,i)')';
    exp = [exp,aveExp];
end

% Density plot
X = data_full(:);
Y = newData(:);

figure (2), clf 
histScatterLog(X,Y,100,10,label('Original', 'BISCUIT-imputed', 'Density plot'),12,'Original_Imputed')
in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

figName = 'Density_Original-Imputed';
save2pdf([figName '.pdf'],gcf,150)
saveas(gcf,figName,'fig')
disp([figName '.pdf + .fig'])
 

% Heatmap 

data1 = data_imputed;
% data1 = importdata('imputed.txt');
% data1 = data1';
zData_row = zscore(data1')';
zData_col = zscore(data1');
% [coefs_row, scores_row, latent_row, tsquared_row, var_exp_row] = pca(zA_hat_row);
% [coefs_col, scores_col, latent_col, tsquared_col, var_exp_col] = pca(zA_hat_col);
pcs = 8;
data_plot = zData_row;

% pcData = zeros(numGenes, numCells);
% pcData_row = scores_row(:,1:pcs);
% pcData_col = scores_col(:,1:pcs);

pcData_row = zData_row;
pcData_col = zData_col;

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
figure(3), clf
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
figure (4), clf 
imagesc(data_plot(inds_ordered_rows,inds_ordered_cols))
% imagesc(pcData_col(inds_ordered_cols,:))
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
    








    
