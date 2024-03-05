function [sigma] = winvar(x, percent, dim)

% Compute the %percent Winsorized variance across dimension dim.
% Implementation according to Wilcox; Robust Estimation
% and Hypothesis Testing. As opposed to matlab's implementation
% the percent is taken on both sides, rather than taking percent/2

if nargin<3
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2
  percent = 0.2;
end

x     = winsorize(x, percent, dim);
n     = size(x, dim);
mu    = mean(x, dim);
sigma = sum(bsxfun(@minus, x, mu).^2, dim)./(n-1);
