function m = bmean(x, dim, k)

% BMEAN computes a skipped-estimator of location based on computing the
% mean after removal of outliers, according to the boxplot rule.
%
% Use as
%   m = bmean(x, dim, k)
%
% Input arguments:
%
%   x   = data matrix
%   dim = scalar, dimension along which to operate, default: first non-singleton
%   k   = scalar threshold for the boxplot rule, default: 1.5

if nargin<3
  k = 1.5;
end

if nargin<2
  dim = find(size(x)>1, 1, 'first');
end

out = boxplotrule(x, dim, k);
x(out) = nan;
m   = nanmean(x, dim);
