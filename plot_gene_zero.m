%% plot genes vs percent of zeros 

addpath(fullfile('~', 'Desktop', 'DATA'));
addpath(fullfile('~','emily_functions'))
load cbcolors_rainbow.mat

% GSM1599494
inputF = importdata('GSM1599494_ES_d0_main.csv');
countmatrix = inputF.data;
cutoff = gene_vs_perZero(countmatrix);
plot(cutoff(:,1),cutoff(:,2), 'LineWidth', 5, 'Color', cbcolors(1,:))
hold on 

% GSM1599495
inputF = importdata('GSM1599495_ES_d0_biorep_techrep1.csv');
countmatrix = inputF.data;
cutoff = gene_vs_perZero(countmatrix);
plot(cutoff(:,1),cutoff(:,2), 'LineWidth', 5, 'Color', cbcolors(2,:))
hold on 

% GSM1599496
inputF = importdata('GSM1599496_ES_d0_biorep_techrep2.csv');
countmatrix = inputF.data;
cutoff = gene_vs_perZero(countmatrix);
plot(cutoff(:,1),cutoff(:,2), 'LineWidth', 5, 'Color', cbcolors(4,:))

grid on
title('Number of genes according to percent of zeros','FontSize',14)
xlabel('Percent of zeros','FontSize',12)
ylabel('Number of genes','FontSize',12)