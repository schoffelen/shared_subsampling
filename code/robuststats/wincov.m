function [sigma] = wincov(x, y, percentx, percenty, dim)

% Compute the %percent Winsorized covariance across dimension dim, between
% the rows or columns of matrices x and y. 
% Implementation according to Wilcox; Robust Estimation
% and Hypothesis Testing. As opposed to matlab's implementation
% the percent is taken on both sides, rather than taking percent/2

if nargin<5
  dim = find(size(x)>1, 1, 'first');
end

if nargin<4 || isempty(percentx)
  percentx = 0.2;
end

if nargin<3 || isempty(percenty)
  percenty = 0.2;
end

% check the input
if ~(size(x,dim)==size(y,dim))
  error('the matrices should have the same number of elements along dimension %d',dim);
end

if any([ndims(x) ndims(y)])>2
  error('the matrices should be at most 2 dimensional');
end

if nargin<2
  percentx = 0.2;
  percenty = 0.2;
end

x     = winsorize(x, percentx, dim);
y     = winsorize(y, percenty, dim);
n     = size(x, dim);
mux   = mean(x, dim);
muy   = mean(y, dim);
% switch dim
%   case 1
%     sigma = bsxfun(@minus, x, mux)'*bsxfun(@minus, y, muy)./(n-1);
%   case 2
%     sigma = bsxfun(@minus, x, mux)*bsxfun(@minus, y, muy)'./(n-1);
%   otherwise
%     error('only value of 1 or 2 is supported for dim');
% end
sigma = sum(bsxfun(@minus, x, mux).*bsxfun(@minus, y, muy), dim)./(n-1);