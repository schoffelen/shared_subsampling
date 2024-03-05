function [omega] = mad(x, dim)

% MAD computes the median absolute deviation across dimension dim.
% Implementation according to Wilcox; Robust Estimation 
% and Hypothesis Testing (page 35)

if nargin<2,
  dim = find(size(x)>1, 1, 'first');
end

pvec  = [dim setdiff(1:ndims(x), dim)];
x     = permute(x, pvec);
x     = demedian(x, 1);
omega = ipermute(nanmedian(abs(x), 1), pvec);
