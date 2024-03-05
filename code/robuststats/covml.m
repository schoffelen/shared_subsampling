function [tmp, mu] = covml(dat, nu, flag)

if nargin<3
  flag = 0;
  % flag = 0: subtract mean
  % flag = 1: scatter only 
end

if nargin<2
  nu = 1;
  % parameter for weighting function
end

[k, nobs] = size(dat);
datorig   = dat;

if ~flag
  mu    = mean(dat,2)*ones(1,nobs);
  error('not yet implemented');
else
  mu    = 0;
end  
tmp      = dat - mu;
sigma    = (tmp*tmp')./nobs;
sigmaold = sigma;

niter = 100;
tmp   = zeros(k,k,niter);
for iter = 1:niter
  s     = sum(conj(dat).*(inv(sigmaold)*dat),1);
  wvec  = (2*k + nu)./(nu + 2.*s);
  sigma = (dat*diag(wvec.*s)*dat')./nobs;
  sigmaold = sigma;
  tmp(:,:,iter) = sigma;
end

