%% plot genes vs percent of zeros 

% add path to Datasets folder

load cbcolors_rainbow.mat

% GSM1599494
inputF = 'GSM1599494_ES_d0_main.csv';
fid = fopen(inputF,'r');
tline = fgetl(fid);
fclose(fid);
columns = cellstr(strvcat(strsplit(tline,'\t')));    % get sample names
totSamps = length(columns); 
% Now that we know # of samples, reopen and get the rest
fid = fopen(inputF,'r');
C = textscan(fid,['%s' repmat('%f',1,totSamps)],'Delimiter','\t','Headerlines',0);
fclose(fid);
genesc = C{1};          % get gene names
countmatrix = [C{2:end}];
cutoff = gene_vs_perZero(countmatrix);
plot(cutoff(:,1),cutoff(:,2), 'LineWidth', 5, 'Color', cbcolors(1,:))
hold on 

% GSM1599495
inputF = 'GSM1599495_ES_d0_biorep_techrep1.csv';
fid = fopen(inputF,'r');
tline = fgetl(fid);
fclose(fid);
columns = cellstr(strvcat(strsplit(tline,'\t')));    % get sample names
totSamps = length(columns); 
% Now that we know # of samples, reopen and get the rest
fid = fopen(inputF,'r');
C = textscan(fid,['%s' repmat('%f',1,totSamps)],'Delimiter','\t','Headerlines',0);
fclose(fid);
genesc = C{1};          % get gene names
countmatrix = [C{2:end}];
cutoff = gene_vs_perZero(countmatrix);
plot(cutoff(:,1),cutoff(:,2), 'LineWidth', 5, 'Color', cbcolors(2,:))
hold on 

% GSM1599496
inputF = 'GSM1599496_ES_d0_biorep_techrep2.csv';
fid = fopen(inputF,'r');
tline = fgetl(fid);
fclose(fid);
columns = cellstr(strvcat(strsplit(tline,'\t')));    % get sample names
totSamps = length(columns); 
% Now that we know # of samples, reopen and get the rest
fid = fopen(inputF,'r');
C = textscan(fid,['%s' repmat('%f',1,totSamps)],'Delimiter','\t','Headerlines',0);
fclose(fid);
genesc = C{1};          % get gene names
countmatrix = [C{2:end}];
cutoff = gene_vs_perZero(countmatrix);
plot(cutoff(:,1),cutoff(:,2), 'LineWidth', 5, 'Color', cbcolors(4,:))

grid on
title('Number of genes according to percent of zeros','FontSize',14)
xlabel('Percent of zeros','FontSize',12)
ylabel('Number of genes','FontSize',12)
