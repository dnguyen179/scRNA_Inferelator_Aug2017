%% random pick 50% GS as priors 

inputF = '/Users/ngu9wm/Desktop/priors/ESCAPE_KC.tsv';
fid = fopen(inputF);
tline=fgetl(fid);
reg_names = cellstr(strvcat(strsplit(tline,'\t')))';
fclose(fid);
totRegs = length(reg_names);

fid = fopen(inputF);
C = textscan(fid,['%s' repmat('%f',1,totRegs)],'Delimiter','\t','Headerlines',1);
fclose(fid);
gene_names = C{1};
totTargs = length(gene_names); 
ncounts = [C{2:end}];
[row, col] = size(ncounts);

prior_matrix = zeros(row, col);
gs_matrix = zeros(row,col);

for i = 1:col
    tf = ncounts(:,i);
    inds = find(tf);
    prior_inds = datasample(inds, round(0.5*length(inds))); 
    ismem = ismember(inds, prior_inds);
    gs_inds = find(ismem==0);
    prior_matrix(prior_inds,i) = ncounts(prior_inds, i);
    gs_matrix(inds(gs_inds),i) = ncounts(inds(gs_inds), i);
end


outFile = fullfile('prior_matrix_50GS.txt');
fout = fopen(outFile,'w');
fprintf(fout,['\t' strjoin(reg_names,'\t') '\n']);
for gene = 1:length(gene_names)
    fprintf(fout,[gene_names{gene} '\t' strjoin(cellstr(num2str(prior_matrix(gene,:)')),'\t') '\n']);
end
fclose(fout);
disp([outFile ' generated.'])
% 
% 
outFile = fullfile('gs_matrix_50GS.txt');
fout = fopen(outFile,'w');
fprintf(fout,['\t' strjoin(reg_names,'\t') '\n']);
for gene = 1:length(gene_names)
    fprintf(fout,[gene_names{gene} '\t' strjoin(cellstr(num2str(gs_matrix(gene,:)')),'\t') '\n']);
end
fclose(fout);
disp([outFile ' generated.'])
