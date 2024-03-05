function [x] = demedian(x, dim)

% Helper function to subtract the median across dimension dim.

if nargin<2,
  dim = find(size(x)>1, 1, 'first');
end

pvec  = [dim setdiff(1:ndims(x), dim)];
x     = permute(x, pvec);
siz   = size(x);
mx    = nanmedian(x, 1);
x     = x - mx(ones(siz(1),1),:,:,:,:,:);
