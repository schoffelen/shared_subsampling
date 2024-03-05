function [b] = bivar(x, kappa, dim)

%Compute biweight midvariance across dimension dim
%Implementation according to Wilcox; Robust Estimation
%and Hypothesis Testing (page 93-96)

if nargin<3,
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2 || isempty(kappa),
  %this is supposed to be a good value for kappa
  kappa = 9;
end

ndim   = ndims(x);
n      = size(x, dim);
repvec = ones(1, ndim);
repvec(dim) = n;

omega = mad(x, dim);
mx    = median(x, dim);
x     = x - repmat(mx, repvec); %de-median x
y     = x./(kappa .* repmat(omega, repvec));
a     = double(abs(y)<1);

num   = sqrt(sum(a.*(x.^2).*((1 - y.^2).^4), dim));
denom = abs(sum(a.*(1 - y.^2).*(1 - 5.*(y.^2)), dim));
b     = (sqrt(n).*(num./denom)).^2;
