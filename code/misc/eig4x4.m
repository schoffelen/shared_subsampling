function [e,d] = eig4x4(in, in2, flag)

if nargin<2
  in2 = false;
end
if nargin<3
  flag = false;
end
if ~isscalar(in2)
  % assume it's a symmetric matrix, to be used for a generalized eigen
  % decomposition
  if nargin<3
    flag = false;
  end
  
  [e,d] = eig4x4(in2, false);
  e     = normoverdim(e,2);
  
  % create the whitening matrix
  d(:,1,1) = 1./sqrt(d(:,1,1));
  d(:,2,2) = 1./sqrt(d(:,2,2));
  d(:,3,3) = 1./sqrt(d(:,3,3));
  d(:,4,4) = 1./sqrt(d(:,4,4));
  
  p     = mtimes4xN(e,d); % use the 4xN version, otherwise the mex file will be used, which operates only correctly on 4x4xN, and not on Nx4x4
  in    = sandwich4x4(permute(p,[1 3 2]), in);
  [e,d] = eig4x4(in, flag);
  e     = normoverdim(mtimes4xN(p,e),2);
  return;
else
  flag = in2;
end

% optimized code to get characteristic polynomial, assuming symmetric
% matrix
c    = charpoly4x4_sym(in);

% get the roots of the polynomial, i.e. the eigenvalues
dtmp = roots_poly4(c, flag);

if flag 
  % only compute the largest
  e = lambda2vec(in, dtmp);
  e = normoverdim(e,2);
  d = dtmp;
else
  
  for k = 1:4
    e(:,:,k) = lambda2vec(in, dtmp(:,k));
    d(:,k,k) = dtmp(:,k);
  end
end

function vec = lambda2vec(in, d)

vec = zeros(size(in,1),4);

in(:,1,1) = in(:,1,1)-d;
in(:,2,2) = in(:,2,2)-d;
in(:,3,3) = in(:,3,3)-d;
in(:,4,4) = in(:,4,4)-d;

% gauss jordan
in(:,1,:)       = in(:,1,:)./in(:,1,1);
in(:,[2 3 4],:) = in(:,[2 3 4],:) - bsxfun(@times,in(:,[2 3 4],1),in(:,1,:));
in(:,2,:)       = in(:,2,:)./in(:,2,2);
in(:,[1 3 4],:) = in(:,[1 3 4],:) - bsxfun(@times,in(:,[1 3 4],2),in(:,2,:));
in(:,3,:)       = in(:,3,:)./in(:,3,3);
in(:,[1 2 4],:) = in(:,[1 2 4],:) - bsxfun(@times,in(:,[1 2 4],3),in(:,3,:));

vec(:,4) = 1;
vec(:,3) = -in(:,3,4);
vec(:,2) = -in(:,2,4);
vec(:,1) = -in(:,1,4);
