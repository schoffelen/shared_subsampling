function [thrlo, thrhi, ix] = percthreshold(x, percent, dim, outflag)

% Helper function to determine the exclusive values of the
% %percent percentile threshold (two-sided) used for 
% trimming, winsorizing etc.
%
% Implementation according to Wilcox; Robust Estimation
% and Hypothesis Testing

if nargin<4 || isempty(outflag)
  outflag = 1;
end

if nargin<3 || isempty(dim)
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2 || isempty(percent)
  percent = 0.2;
end

n    = size(x, dim);
g    = floor(percent*n);

pvec  = [dim setdiff(1:ndims(x), dim)];
[sx, ix] = sort(permute(x, pvec), 1);
thrlo = ipermute(sx(g+1,:,:,:,:,:,:,:),pvec);
thrhi = ipermute(sx(n-g,:,:,:,:,:,:,:),pvec);

switch outflag
  case 0
    ix    = [];
  case 1
    % indices
    ix    = ipermute(ix((g+1):(n-g),:,:,:,:,:,:,:),pvec);
  case 2
    % booleans
    [~, ix] = sort(ix, 1);
    ix      = ipermute(ix>g&ix<=n-g,pvec);
  otherwise
end
