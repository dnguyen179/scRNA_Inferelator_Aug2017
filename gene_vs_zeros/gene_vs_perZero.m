function[numgene_perzeros] = gene_vs_perZero(inputData)
%% Plot percent genes vs sparsity of gene row 

countmatrix = inputData;
[Nrow, Ncol] = size(countmatrix);

perZeros = zeros(Nrow,1);
count = 0; 
for i = 1:Nrow
    for j = 1:Ncol
        if countmatrix(i,j) == 0
            count = count + 1;
        end
    end     
    perZeros(i,1) = (count/Ncol)*100;
    count = 0; 
      
end

cutoff = [100,2];
numGenes = 0;
for i = 1:100
    cutoff(i,1) = i;
    for j = 1:Nrow
        if perZeros(j,1) >= i
            numGenes = numGenes + 1;
            cutoff(i, 2) = numGenes;
        end 
    end
    numGenes = 0;
end

numgene_perzeros = cutoff;
