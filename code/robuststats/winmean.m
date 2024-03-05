function [mu] = winmean(x, percent, dim)

%Compute the %percent Winsorized mean across dimension dim.
%Implementation according to Wilcox; Robust Estimation
%and Hypothesis Testing

if nargin<3,
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2,
  percent = 0.2;
end

x  = winsorize(x, percent, dim);
mu = mean(x, dim);
