function [mat_inf_new] = get_kmeans_clustering(...
    data,...             % COMPLETE data matrix
    x_labelsc,...        % cell string for X labels
    y_labelsc,...        % cell string for Y labels
    dataset_name,...     % string for dataset name
    x_clusterings,...    % if unknown, set to [].  
    ...                  % if specified, x_clustering denotes cluster membership
    desired_x_clusters,...  % number of x clusters desired.  
    ...                  % if clusters already specified, set to 0 or ''
    y_clusterings,...    % if unknown, set to [].  if specified, y_clustering
    ...                  % denotes cluster membership
    desired_y_clusters,...  % number of y clusters desired.  
    ...                  % if y clusters already specified, set to 0 or ''
    data_scaling,...     % choose among the following: 'Zscore', 'log2',
    ...                  % or 'FC' (FC = fold-change), 'none' (no scaling)
    distance_metric,...       % distance used to calculate similarity matrix
    ...                  % choose 'Correlation' or 'Euclidean';
    mat_inf,...          % Name Base for Outputs (should be a string)
    general_function_path,...   % path to general functions / if not in same folder, set to ''
    clustering_function_path)  % path to clustering functions / if in same folder, set to ''
% %% function [mat_inf] = get_clustering_solutions_function(...
% %     data,...             % COMPLETE data matrix
% %     x_labelsc,...        % cell string for X labels
% %     y_labelsc,...        % cell string for Y labels
% %     dataset_name,...     % string for dataset name
% %     x_clusterings,...    % if unknown, set to [].  
% %     ...                  % if specified, x_clustering denotes cluster membership
% %     desired_x_clusters,...  % number of x clusters desired.  
% %     ...                  % if clusters already specified, set to 0 or ''
% %     y_clusterings,...    % if unknown, set to [].  if specified, y_clustering
% %     ...                  % denotes cluster membership
% %     desired_y_clusters,...  % number of y clusters desired.  
% %     ...                  % if y clusters already specified, set to 0 or ''
% %     data_scaling,...     % choose among the following: 'Zscore', 'none' (no scaling)
% %     distance_metric,...       % distance used to calculate similarity matrix
% %     ...                  % choose 'Correlation' or 'Euclidean';
% %     mat_inf,...          % Name Base for Outputs (should be a string)
% %     general_function_path,...   % path to general functions / if not in same folder, set to ''
% %     clustering_function_path)  % path to clustering functions / if in same folder, set to ''
% % Function to cluster a COMPLETE data matrix with kmeans

% cd ~/Mick/Matlab
% load GSEA/MDDC_GSEA_biocarta_1000.mat
% general_function_path = path_to_emsFXNs;   % path to general functions / if not in same folder, set to ''
    
%% parameters for debugging, if necessary
% % maximum fraction of missing data in 
% % columns' direction
% max_fraction_missing_x_direction = .2;
% % rows' direction
% max_fraction_missing_y_direction = 0;
% % x-direction clustering
% x_clusterings = currcond_nums;  % do not want to cluster in x direction
% desired_x_clusters = [];
% % y-direction clustering
% y_clusterings = [];             % to be determined
% desired_y_clusters = [4];     
% % Data Scaling - 'Zscore', 'log2', or 'FC' (FC = fold-change), 'none' (no
% % scaling)
% data_scaling = 'log2';
% % Distance Metric - Euclidean or Correlation
% distance_metric = 'Correlation';
% % p-value option, if 1 --> p-val incorporated into distance metric
% pval_opt = 1;                   % if 0 --> no p-val incorporation
% % Name Base for Output
% mat_inf = 'ClusteringExample';
% 
% % Set Dataset Variables
% x_labelsc = cellstr(num2str(currlivernums));   % cell string
% y_labelsc = cellstr(strvcat(currlabels));                         % cell string
% data = currdata;                                % 0 or -9 denote missing data
% dataset_std = currdata_std;
% dataset_name = 'Basal pY';              
%
% % Set Path 
% % path to my general functions
% general_function_path = 'C:\Documents and Settings\emiraldi\My Documents\MATLAB\emily_functions';
% % path to my clustering functions
% clustering_function_path = 'C:\Documents and Settings\emiraldi\My Documents\MATLAB\emily_functions\clustering';

%% Get paths to functions, if they're not in the same file
if general_function_path
    addpath(general_function_path)
end
if clustering_function_path
    addpath(clustering_function_path)
end

%% Apply data transformations
if ismember({'Z-score'},data_scaling)  % for z-score scaling
    plotdata = data';
    plotdata = zscore(plotdata);
elseif ismember({'sqrt(Log2)'},data_scaling)
    error('Not implemented')
elseif ismember({'Log2'},data_scaling)
    error('Not implemented')
else  % no scaling fpkm
    plotdata = data';
end

rodata = plotdata';
codata = plotdata;
plotdata = plotdata';
[rownum colnum] = size(plotdata);

%% Clustering
% X / column direction: 
% run Kmeans algorithm, if necessary
if desired_x_clusters
    kids_cols = kmeans(codata,desired_x_clusters,'Distance',distance_metric);
else 
    kids_cols = ones(colnum,1);
end
% cluster within clusters, using hierarchical clustering
inds_ordered = zeros(colnum,1);
ex_ids = unique(kids_cols);
clust_num = length(ex_ids);
start_spots = zeros(clust_num,1);
clust_sizes = zeros(clust_num,1);
start = 1;
for clust = 1:clust_num          
    clust_ind = find(kids_cols == ex_ids(clust));
    if length(clust_ind) > 1
        pdis = pdist(codata(clust_ind,:));
        link = linkage(pdis,'average');
        [h t horder] = dendrogram(link,0);
    else
        horder = 1;
    end
    clust_size = length(clust_ind);
    clust_sizes(clust) = clust_size;
    % order indices
    inds_ordered(start:start+clust_size-1) = clust_ind(horder);
    start = clust_size+start;
    start_spots(clust) = start;
end    

% Y / row direction
% run Kmeans algorithm, if necessary
if desired_y_clusters
    [kids_rows centers sumd,D] = kmeans(rodata,desired_y_clusters,'Distance',distance_metric,'MaxIter',10000);
else
    kids_rows = ones(rownum,1);
    centers = mean(rodata);
end
% cluster within clusters
inds_ordered_rows = zeros(rownum,1);
ex_ids_rows = unique(kids_rows);
clust_num_rows = length(ex_ids_rows);
start_spots_rows = zeros(clust_num_rows,1);
clust_sizes_rows = zeros(clust_num_rows,1);
start_rows = 1;
for clust = 1:clust_num_rows
    clust_ind_rows = find(kids_rows == ex_ids_rows(clust));
    if length(clust_ind_rows) > 1
        pdis = pdist(rodata(clust_ind_rows,:));
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

mat_inf_new = [mat_inf '_' distance_metric '_' data_scaling  '_XClust' ...
    num2str(clust_num) '_YClust' num2str(clust_num_rows)];

clustering_type = 'K-means';

save([mat_inf_new '.mat'],...
    'D',...
    'centers',...                            4x24                  768  double              
    'clust',...                              1x1                     8  double              
    'clust_num',...                          1x1                     8  double              
    'clust_num_rows',...                     1x1                     8  double              
    'clust_size',...                         1x1                     8  double              
    'clust_size_rows',...                    1x1                     8  double              
    'clust_sizes',...                        2x1                    16  double              
    'clust_sizes_rows',...                   4x1                    32  double              
    'clustering_function_path',...           1x35                   70  char 
    'clustering_type',...
    'codata',...                            24x292               56064  double              
    'colnum',...                             1x1                     8  double              
    'data',...                             292x24                56064  double              
    'data_scaling',...                       1x4                     8  char                
    'dataset_name',...                       1x4                     8  char                
    'desired_x_clusters',...                 1x1                     8  double              
    'desired_y_clusters',...                 1x1                     8  double              
    'distance_metric',...                    1x11                   22  char                
    'ex_ids',...                             2x1                    16  double              
    'ex_ids_rows',...                        4x1                    32  double              
    'general_function_path',...              1x24                   48  char                
    'inds_ordered',...                      24x1                   192  double              
    'inds_ordered_rows',...                292x1                  2336  double              
    'kids_cols',...                         24x1                   192  double              
    'kids_rows',...                        292x1                  2336  double              
    'mat_inf',...
    'mat_inf_new',...1x63                  126  char                
...    pdis                               1x1770              14160  double              
    'plotdata',...                         292x24                56064  double              
    'rodata',...                           292x24                56064  double              
    'rownum',...                             1x1                     8  double              
    'start_rows',...                         1x1                     8  double              
    'start_spots',...                        2x1                    16  double              
    'start_spots_rows',...                   4x1                    32  double              
    'x_clusterings',...                      0x0                     0  double              
    'x_labelsc',...                         24x1                  3132  cell                
    'y_clusterings',...                      0x0                     0  double              
    'y_labelsc');%,...                        292x1                 35782  cell                
disp(['Generated: ' mat_inf_new  '.mat'])