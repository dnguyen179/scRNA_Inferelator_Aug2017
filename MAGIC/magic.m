%% MAGIC normalization 

addpath(fullfile('~','Miraldi_functions'))

% Load data and gene/cell names - input matrix should be cells x genes 
inputF = importdata('GSM1599497_ES_d2_LIFminus.csv');
OriginalData = inputF.data;  % genesxcells matrix

data = OriginalData'; % cellsxgenes matrix as input for MAGIC

% run MAGIC
npca = 20; % ususally between 10 and 200
ka = 5; % can be smaller, eg 3 
k = 15; % can be smaller, eg 9 (3*ka)
t = 6; % usually between 6 and 12, smaller ka/k requitres bigger t
lib_size_norm = true; % library size normalize
log_transform = false; % log transform, some data requires this

% remember to change the values of epsilon, rescale_to, pseudo_count
data_imputed = run_magic(data, t, 'npca', npca, 'ka', ka, 'k', k, ...
    'lib_size_norm', lib_size_norm, 'log_transform', log_transform);

newData = transpose(data_imputed); %genesxcells matrix for subsequent analysis

%% Rescale newData

[numGenes, numCells] = size(newData);

libsize  = sum(transpose(OriginalData),2);
OriginalDataNorm = bsxfun(@rdivide, transpose(OriginalData), libsize) * median(libsize);
OriginalDataNorm = OriginalDataNorm';

nonZeroOriData = nonzeros(OriginalDataNorm);
nonZeroNewData = nonzeros(newData);
med_oriData = median(nonZeroOriData);
med_newData = median(nonZeroNewData);
med_dev = med_oriData - med_newData;

mad_oriData = mad(nonZeroOriData);
mad_newData = mad(nonZeroNewData);
rescale = mad_oriData/mad_newData;

Scaled_Data = zeros(numGenes, numCells);
for i = 1:numGenes
    for j = 1:numCells
        Scaled_Data(i,j) = (newData(i,j) + med_dev)*rescale;
    end
end

save('gsm97_magic.mat', 'Scaled_Data');

%% Imputed vs original gene expression 

libsize  = sum(transpose(OriginalData),2);
OriginalDataNorm = bsxfun(@rdivide, transpose(OriginalData), libsize) * median(libsize);
OriginalDataNorm = OriginalDataNorm';

X = log(OriginalDataNorm(:)+1);
Y= log(Scaled_Data(:)+1);

% Distribution histogram
figure (1), clf
fig1 = subplot(2,1,1);
hist(fig1,X,75);
axis tight
title('Distribution of log gene expression of normalized raw data','FontSize',14)
xlabel('Log of gene expression','FontSize',12)
ylabel('Abundance','FontSize',12)

fig2 = subplot(2,1,2);
hist(fig2,Y,75);
axis tight
title('Distribution of log gene expression of MAGIC imputed data','FontSize',14)
xlabel('Log of gene expression','FontSize',12)
ylabel('Abundance','FontSize',12)

in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

figName = 'Dist_Original-Imputed';
save2pdf([figName '.pdf'],gcf,150)
saveas(gcf,figName,'fig')
disp([figName '.pdf + .fig'])

% Density plot
figure (2), clf 
histScatterLog(X,Y,100,10,label('Original', 'MAGIC', 'MAGIC-Density plot'),16, 'Original_Imputed')
in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');

figName = 'Density_Original-Imputed';
save2pdf([figName '.pdf'],gcf,200)
saveas(gcf,figName,'fig')
disp([figName '.pdf + .fig'])
 
%% Generate heatmap

%z-score  
zScaled_Data_row = zscore(Scaled_Data')';
zScaled_Data_col = zscore(Scaled_Data');

pcs = 4; % number of desired clusters 
pcData = zScaled_Data_row;
pcData_row = zScaled_Data_row;
pcData_col = zScaled_Data_col; 

% k-means
kids_rows = kmeans(pcData_row,pcs);
kids_cols = kmeans(pcData_col,pcs);

% cluster within clusters - row direction
inds_ordered_rows = zeros(numCells,1);
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
imagesc(pcData(inds_ordered_rows,inds_ordered_cols))
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
in=input('Adjust figure dimensions on your screen, and, once satisfied, press any key to resume.');
figName = 'HeatmapClusterExp_rowcol';
save2pdf([figName '.pdf'],gcf,150)
saveas(gcf,figName,'fig')
disp([figName '.pdf + .fig'])





