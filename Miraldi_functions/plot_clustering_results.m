function plot_clustering_results(...
    mat_inf,...     % specify the .mat file with clustering results
    output_folder,...   % output folder for figures and such
    ...             % Example, if results in cluster.mat, type 'cluster'
    axis_fontsize,...           % font size of axis
    title_fontsize,...          % font size of title
    papersizeX,...              % size of paper in x-direction (inches) (for pdf file)
    papersizeY,...              % size of paper in y-direction (inches)
    degree_x_label_rotation,... % degree to rotate x labels
    dot_size,...                % size of dot which will denote missing data 
    general_function_path)   % path to general functions / if not in same folder, set to ''
% function to plot clustering results stored in .mat file created by
% get_clustering_solutions.m

%% parameters for debugging if necessary
% mat_inf = 'ClusteringExample_log2_pval_XClust4_YClust4_Xmiss20_Ymiss0';
% 
% [rows cols] = size(plotdata);
% 
% % Figure parameters
% axis_fontsize = 10;             % font size of axis
% title_fontsize = 16;            % font size of title
% papersizeX = 11;                % size of paper in x-direction (inches)
% papersizeY = 11;               % size of paper in y-direction (inches)
% degree_x_label_rotation = 90;    % degree to rotate x labels
% missing_data_option = 'dot';    % missing data represented as a 'dot' or 'square'
% % --> this option doesn't work with my version of Matlab, 2011b or higher
% % needed
% dot_size = 5;

%% Load Data
load([mat_inf '.mat'])
[rows cols] = size(plotdata);

dist_label = 'Squared Euclidean';
data_scaling = 'Log2';

% mat_inf = fullfile(output_folder,mat_inf);

%% Get paths to functions, if they're not in the same file
if general_function_path
    addpath(general_function_path)
end

% find out whether data scale should be zero or 1 centered
rodatavector = rodata(:);
min_data = min(rodatavector);

%% DATA Figure -- Limit Color Max
figure
imagesc(plotdata(inds_ordered_rows,inds_ordered))
set(gca,'YTick',1:rows,'YTickLabel',strvcat(y_labelsc{inds_ordered_rows}))
xticklabel_rotate([1:cols],degree_x_label_rotation,cellstr(strvcat(x_labelsc{inds_ordered})),...
    'interpreter','none','FontSize',axis_fontsize,'FontWeight','Bold')

set(gca,'FontSize',axis_fontsize,'LineWidth',2,'FontWeight','Bold')
ax = axis();
hold on
for clust = 1:clust_num_rows
    plot([ax(1) ax(2)], [start_spots_rows(clust) start_spots_rows(clust)]-.5,'k',...
        'LineWidth',2)
end
hold on
for clust = 1:clust_num
    plot([start_spots(clust) start_spots(clust)]-.5,[ax(3) ax(4)],'k',...
        'LineWidth',2)
end

% if length(find(ismember(missing_data_option,'dot'))) > 0
%     [ro co] = find(data(inds_ordered_rows,inds_ordered) == -9);  % find missing data for plotting
%     plot(co,ro,'.','Color',[.75 .75 .75],'MarkerSize',dot_size)
% % else
% %     mask = make_mask(data(inds_ordered_rows,inds_ordered));
% %     set(gcf,'AlphaData',mask)
% %     set(gca,'YTickLabel',[],'XTickLabel',[],'PlotBoxAspectRatio',[2 4 1],'YTick',0,'XTick',0)
% %     whitebg([0.4773    0.5710    0.6023])
% end

title(['Kmeans (' dist_label ', ' data_scaling, ')'],...
    'Fontsize',title_fontsize)
if min_data >=0
    colormap redblue
    set(gca,'CLim',[0 2])
    colorbar('Location','EastOutside','YTick',[0:2],...
        'YTickLabel',strvcat('0','1','>2'),'FontSize',axis_fontsize,...
        'FontWeight','Bold','LineWidth',2)        
else
    colormap redblue
    set(gca,'CLim',[-2 2])
    colorbar('Location','EastOutside','YTick',[-2:1:2],'YTickLabel',...
        strvcat('<-2',num2str([-1:1]'),'>2'),'FontSize',axis_fontsize,...
        'FontWeight','Bold','LineWidth',2)    
end

fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [papersizeX papersizeY]);
% print('-painters','-dpdf','-r1800',...
%     [mat_inf '_sat.pdf'])
saveas(gcf,[mat_inf '_sat'],'tiff')
saveas(gcf,[mat_inf '_sat'],'fig')

%% DATA Figure -- unsaturated colors
figure
imagesc(plotdata(inds_ordered_rows,inds_ordered))
set(gca,'YTick',1:rows,'YTickLabel',strvcat(y_labelsc{inds_ordered_rows}))
xticklabel_rotate([1:cols],degree_x_label_rotation,cellstr(strvcat(x_labelsc{inds_ordered})),...
    'interpreter','none','FontSize',axis_fontsize,'FontWeight','Bold')

colormap redblue
set(gca,'FontSize',axis_fontsize,'LineWidth',2,'FontWeight','Bold')
ax = axis();
hold on
for clust = 1:clust_num_rows
    plot([ax(1) ax(2)], [start_spots_rows(clust) start_spots_rows(clust)]-.5,'k',...
        'LineWidth',2)
end
hold on
for clust = 1:clust_num
    plot([start_spots(clust) start_spots(clust)]-.5,[ax(3) ax(4)],'k',...
        'LineWidth',2)
end

% denote missing data
% if length(find(ismember(missing_data_option,'dot'))) > 0
%     plot(co,ro,'.','Color',[.75 .75 .75],'MarkerSize',dot_size)
% % else
% %     mask = make_mask(data(inds_ordered_rows,inds_ordered));
% %     set(gcf,'AlphaData',mask)
% %     set(gca,'YTickLabel',[],'XTickLabel',[],'PlotBoxAspectRatio',[2 4 1],'YTick',0,'XTick',0)
% %     whitebg([0.4773    0.5710    0.6023])
% end

title(['K-means (' dist_label ', ' data_scaling ')'],'Fontsize',title_fontsize)
    %...', BIC = ', num2str(rowund(bic)), ', ' dist_label])       
if min_data >=0
    colormap jet
    colorbar('Location','EastOutside','FontSize',axis_fontsize,...
        'FontWeight','Bold','LineWidth',2)    
else
    colormap redblue
    maxabs = max(max(max(plotdata)),-min(min(plotdata)));    
    set(gca,'CLim',[-maxabs maxabs],'FontSize',axis_fontsize)
    colorbar('Location','EastOutside',...
        'FontSize',axis_fontsize,'FontWeight','Bold','LineWidth',2)%,'YTick',[-2:1:2],'YTickLabel',...
end
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [papersizeX papersizeY]);
% print('-painters','-dpdf','-r1800',...
%     [mat_inf  '.pdf'])
saveas(gcf,[mat_inf ],'tiff')
saveas(gcf,[mat_inf  ],'fig')

% %% rows sims figure
% figure
% imagesc(Sim_rows(inds_ordered_rows,inds_ordered_rows))
% set(gca,'YTick',1:rows,'YTickLabel',strvcat(y_labelsc{inds_ordered_rows}))
% set(gca,'FontSize',axis_fontsize,'FontWeight','Bold')
% xticklabel_rotate([1:rows],degree_x_label_rotation,cellstr(strvcat(y_labelsc{inds_ordered_rows})),...
%     'interpreter','none','FontSize',axis_fontsize,'FontWeight','Bold')
% %set(gca,'XTick',find(mcs_all_bol(inds_ordered_rows)),'XTickLabel','*')
% % xlabel(dist_label)
% set(gca,'FontSize',axis_fontsize,'LineWidth',2,'FontWeight','Bold')
% colorbar('Location','EastOutside','FontSize',axis_fontsize,'FontWeight','Bold',...
%     'LineWidth',2)
% 
% hold on
% ax = axis();
% for clust = 1:clust_num_rows
%     plot([ax(1) ax(2)], [start_spots_rows(clust) start_spots_rows(clust)]-.5,'k',...
%         'LineWidth',2)
% end
% for clust = 1:clust_num_rows
%     plot([start_spots_rows(clust) start_spots_rows(clust)]-.5,[ax(3) ax(4)],'k',...
%         'LineWidth',2)
% end
% title(['K-means (' dist_label ')'],'Fontsize',title_fontsize)
% colormap hot
% 
% 
% if length(dist_label) == length('Correlation')
%     set(gca,'CLim',[-1 1])    
% end
% 
% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [papersizeX papersizeY]);
% print('-painters','-dpdf','-r1800',...
%     [mat_inf  '_row.pdf'])
% saveas(gcf,[mat_inf '_row'],'tiff')
% saveas(gcf,[mat_inf  '_row'],'fig')
% 
% %% VARS sims figure
% figure
% imagesc(Sim_cols(inds_ordered,inds_ordered))
% set(gca,'YTick',1:cols,'YTickLabel',strvcat(x_labelsc{inds_ordered}))
% set(gca,'FontSize',axis_fontsize,'LineWidth',2,'FontWeight','Bold')
% xticklabel_rotate([1:cols],degree_x_label_rotation,cellstr(strvcat(x_labelsc{inds_ordered})),...
%     'interpreter','none','FontSize',axis_fontsize,'FontWeight','Bold')
% colormap redblue
% colorbar('Location','EastOutside','FontSize',axis_fontsize,'FontWeight','Bold',...
%     'LineWidth',2)
% colormap hot
% hold on
% ax = axis();
% for clust = 1:clust_num
%     plot([ax(1) ax(2)], [start_spots(clust) start_spots(clust)]-.5,'k',...
%         'LineWidth',2)
% end
% for clust = 1:clust_num
%     plot([start_spots(clust) start_spots(clust)]-.5,[ax(3) ax(4)],'k',...
%         'LineWidth',2)
% end
% 
% if length(dist_label) == length('Correlation')
%     set(gca,'CLim',[-1 1])    
% end
% 
% title(['K-means (' dist_label ')'],'Fontsize',title_fontsize)
% 
% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [papersizeX papersizeY]);
% print('-painters','-dpdf','-r1800',...
%     [mat_inf  '_col.pdf'])
% saveas(gcf,[mat_inf '_col'],'tiff')
% saveas(gcf,[mat_inf '_col'],'fig')
