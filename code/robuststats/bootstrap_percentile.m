function [lower,upper] = bootstrap_percentile(x, varargin)

dim   = ft_getopt(varargin, 'dim', []);
nboot = ft_getopt(varargin, 'nboot', 1000);
alpha = ft_getopt(varargin, 'alpha', 0.05);
est   = ft_getopt(varargin, 'est',   'trimmean');

if isempty(dim)
  dim = find(size(x)>1,1,'first');
end

if ndims(x)>2
  error('at the moment a dimensionality of the data > 2 is not supported');
end

if dim~=2
  error('at the moment only a dim of 2 is supported');
end
nobs = size(x,2);

lb = round(alpha*nboot/2);

L(size(x,1),lb+1) = nan;
U(size(x,1),lb+1) = nan;

switch est
  case 'trimmean'
    % additional parameters
    percent = ft_getopt(varargin, 'percent', 0.2);
  otherwise
    error('only trimmean as an estimator is currently supported');
end

for k = 1:nboot
  sel = ceil(rand(1,nobs));
  
  switch est
    case 'trimmean'
      estb = trimmean(x(:,sel),percent,2);
      
  end
  if k<=lb
    % store the individual sample, because pruning can only take place,
    % once we have more than alph*nboot replicates
    L(:,k) = estb;
    U(:,k) = estb;
  else
    L = sort([L(:,1:lb) estb],2,'ascend');
    U = sort([U(:,1:lb) estb],2,'descend');
  end 
end

upper = U(:,lb);
lower = L(:,lb);