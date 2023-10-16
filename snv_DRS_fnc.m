function [snv_X] = snv_DRS(X)
% A standard normal variate (SNV) normalization was used to process spectral data.
% Using SNV normalization, the mean of each spectrum was set to zero and the 
% standard deviation was set to one.

% Input: X: (m x n); spectral data where m = number of observations and n = number of features.

% Output: snv_X: (m x n); snv transformed data.

% Celina L. Li, Sept 2021.

[m, n] = size(X);

snv_X = ( X - mean(X,2).*ones(1,n) ) ./ ( std(X,0,2).*ones(1,n) );

% end