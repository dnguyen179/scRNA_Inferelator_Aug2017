function [normalized_data] = libsize_norm(input_data)

%% Library size normalization 
DataMatrix = input_data;
libsize  = sum(transpose(DataMatrix),2);
DataNorm = bsxfun(@rdivide, transpose(DataMatrix), libsize) * 10^6;
DataNorm = DataNorm';

normalized_data = DataNorm;
