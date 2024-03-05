function out = madmedianrule(x, dim, k)

% MADMEDIANRULE identifies outliers according to the madmedianrule described in
% Wilcox 3.13.4. 
%
% Use as:
%   out = madmedianrule(x, dim, k)
%
% Input arguments:
%   x   = data matrix
%   dim = dimension along which to operate, default: first non-singleton dimension
%   k   = threshold value, default: 2.24
%
% Output arguments:
%   out = binary mask, size(out) = size(x), true for outliers
% 
% Nans are treated as missing values

if nargin<3
  k = 2.24;
end

if nargin<2
  dim = find(size(x)>1, 1, 'first');
end

m      = nanmedian(x, dim);
repvec = ones(1,ndims(x));
repvec(dim) = size(x,dim);
out    = abs(x-repmat(m,repvec))./(repmat(mad(x,dim),repvec)./0.6745)>k;
