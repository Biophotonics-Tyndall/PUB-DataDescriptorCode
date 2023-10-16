function [y_norm] = AUC_norm(x,y)
% To perform area-under-curve normalization, so all AUCs are 1.
% Utilize Matlab trapz function

% Input: x: (1 x n); scalar, spectral wavelength
%        y: (m x n); spectral data where m = number of observations and n = number of features.

% Output: y_norm: (m x 1); AUCs.

% Celina L. Li, Sept 2021.

area = trapz(x,y,2);

y_norm = y./area;

% end
