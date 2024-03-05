function [beta, R, r] = weighted_regression(y,x,w,wflag)

% Perform weighted_regression of y onto x using weights w
% Outputs: beta=beta weights, R=r-square value, r=residuals
% Use as [beta, R, r] = weighted_regression(y,x,w)

if nargin<4,
  wflag = 0;
end

if ~wflag,
  %compute output using voxel specific weights
  [m,n] = size(w);
  if m~=n && (m==1 || n==1)
    w = repmat(w(:), [1 size(y,2)]);
  end
  wy    = w.*y; %weight the data
  if size(x,2)==1,
    num   = x'*wy;
    denom = sum(repmat(abs(x).^2, [1 size(y,2)]).*w);
    beta  = num./denom;
    r     = y - x*beta;
  elseif size(x)==size(y)
    num   = sum(conj(x).*wy);
    denom = sum(abs(x).^2.*w);
    beta  = num./denom;
    r     = y - x*diag(beta);
  else
    num   = x'*wy; %Nregressors x Nvoxels
    denom = transpose(abs(x).^2)*w;
    beta  = num./denom;
    r     = y - (x*beta)./size(x,2);
  end
  clear num denom;
  
  num   = sum(abs(w.*r).^2); %weighted residuals
  denom = sum(abs(wy).^2);
  R     = 1 - num./denom;
else
  %compute output using single set of weights

  [m,n] = size(w);
  if ~any([m n]==1) && m~=n, 
    error('w should be a vector or a diagonal matrix');
    % FIXME make check more robust
  else
    w = diag(w);
  end

  wy   = w*y;
  beta = (x'*w*x)\(x'*wy);
  r    = y - x*beta; 
  
  num   = sum(abs(w*r).^2);
  denom = sum(abs(wy).^2);
  R     = 1 - num./denom;
end

