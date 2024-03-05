function [q1, q2] = idealf(x, dim)

% IDEALF computes the 'ideal fourths' according to Wilcox paragraph 3.12.5,
% accounting for nans
%
% Use as
%   [q1, q2] = idealf(x, dim)
%
% Input arguments:
%   x   = matrix containing the data
%   dim = dimension along which to operate, if not specified the first
%         non-singleton dimension is used.
%
% Nans in the input are treated as missing data.

if nargin<2
  dim = find(size(x)>1, 1, 'first');
end

% reshape the data into 2D
siz  = size(x);
ix   = setdiff(1:ndims(x),dim);
x    = permute(x, [ix dim]);
x    = reshape(x, [prod(siz(ix)) siz(dim)]);

% pre compute some variables
s = sort(x,2);
n = sum(isfinite(x),2);
j = floor((n/4)+5/12);
h = n/4+5/12-j;
k = n-j+1;

% loop over the rows
q1 = nan(size(x,1),1);
q2 = nan(size(x,1),1);
for i = 1:size(x,1)
  if n(i)>0
    q1(i) = (1-h(i))*s(i,j(i)) + h(i)*s(i,j(i)+1);
    q2(i) = (1-h(i))*s(i,k(i)) + h(i)*s(i,k(i)-1);
  end
end

% reshape back
q1 = ipermute(reshape(q1, [siz(ix) 1]), [ix dim]);
q2 = ipermute(reshape(q2, [siz(ix) 1]), [ix dim]);
