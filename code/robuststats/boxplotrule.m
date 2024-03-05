function out = boxplotrule(x, dim, k)

% BOXPLOTRULE identifies outliers according to the boxplotrule described in
% Wilcox 3.13.2. It uses the ideal fourths, rather than the interquartile
% range from the actual data.
%
% Use as:
%   out = boxplotrule(x, dim)
%
% Input arguments:
%   x   = data matrix
%   dim = dimension along which to operate
%
% Output arguments:
%   out = binary mask, size(out) = size(x), true for outliers
% 
% Nans are treated as missing values

if nargin<3
  k = 1.5;
end

if nargin<2
  dim = find(size(x)>1, 1, 'first');
end

[q1, q2] = idealf(x, dim);

dq     = q2-q1;
repvec = ones(1,ndims(x));
repvec(dim) = size(x,dim);
out    = x<repmat(q1-k.*dq,repvec)|x>repmat(q2+k.*dq,repvec);
