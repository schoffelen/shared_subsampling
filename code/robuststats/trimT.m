function [T] = trimT(x, percent, dim)

if nargin<3,
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2,
  percent = 0.2;
end

m = trimmean(x, percent, dim);
s = trimsem(x, percent, dim);
T = m./s;
