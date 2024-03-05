function [mx] = wpow(x, dim)

%Compute robust estimate of power based on a iterative weighted procedure,
%according to Chave et al. 1987

if nargin<2,
  dim = find(size(x)>1, 1, 'first');
end

if ndims(x)>2, error('dimensionality > 2 is not supported'); end

pvec  = [dim setdiff(1:ndims(x), dim)];
x     = permute(x, pvec);
n     = size(x, 1);
m     = median(x, 1);
r     = x - m(ones(n,1),:);
d     = mad(x, dim)./0.5; %scaling factor FIXME
b     = 2; %beta in (27) in paper FIXME
w     = thomsonweight(r./d, b);

for k = 2:100
  wold(:,k-1) = w;
  mold(1,k-1) = m; 
  m         = diag(1./(sum(w.^2))) .* w' * r;
  r         = x - m(ones(n,1),:);
  w         = thomsonweight(r./d, b);
end

keyboard

function [w] = thomsonweight(r, b);

%Compute weighting function according to reference formula 27,
%but not taking the absolute for x

w = exp(-exp(b.*(abs(r)-b)));
