function [coh, coh0, w, varargout] = wcoh(y, x, dim, w)

%[coh,coh0,w] = wcoh(y,x,dim,w);
%Compute robust estimate of coherence based on a iterative weighted procedure,
%according to Chave et al. 1987


if nargin<3 || isempty(dim),
  dim = find(size(x)>1, 1, 'first');
end

if nargin<4,
  w   = ones(size(y)); %weights 
end

if ndims(x)>2,
  error('dimensionality > 2 is not supported'); 
end

pvec  = [dim setdiff(1:ndims(x), dim)];
x     = permute(x, pvec);
y     = permute(y, pvec);
w     = permute(w, pvec);

sizx  = size(x);
sizy  = size(y);
rvecx = [sizx(1) prod(sizx(2:end))];
rvecy = [sizy(1) prod(sizy(2:end))];
x     = reshape(x, rvecx);
y     = reshape(y, rvecy); %reshape into Nobs x Nvox
w     = reshape(w, rvecy);

n     = sizx(1);            %number of observations
%b     = sqrt(2.*log(2.*n)); %parameter for weighting function (eq 33)
b     = sqrt(2.*log(n./3)); %parameter for weighting function corresponding with 97.5th quantile
%b     = sqrt(2.*log(n./5.5));[;
%disp(b)

numiter = 1000;
thresh  = 1e-15;
ok      = logical(zeros(1,sizy(2)));
okold   = logical(zeros(1,sizy(2)));
Rold    = zeros(1,sizy(2)) + inf;

for k = 1:numiter
  %[beta, R, r] = weighted_regression(y(:,~ok), x, w(:,~ok));
  [beta, R, r] = weighted_regression(y(:,~ok), x(:,~ok), w(:,~ok));
  
  d        = mad(abs(r), 1)./0.44845; %scaling parameter
  w(:,~ok) = thomsonweight(r./d(ones(1,n),:), b);

  ok(~ok)   = abs(R-Rold(~ok)) < thresh;
  
  Rold(~okold)   = R;
  okold          = ok;

  if sum(ok) == sizy(2),
    break
  end

  if k==1,
    coh0  = R;
    %beta0 = beta;
  end
end
coh = Rold;

if nargout>3,
  varargout{1} = cohZtransform(coh, sum(w.^2));
  varargout{2} = cohZtransform(coh0, n);
else
  varargout = {};
end

%----------------------------------
function [w] = thomsonweight(r, b)

%Compute weighting function according to reference formula 27,

w = exp(-exp(b.*(abs(r)-b)));

%------------------------------------
function [Z] = cohZtransform(coh, n)

%Ztransform according to Jarvis and Mitra
%note that coh is already squared and absolute

beta = 23./20;
q    = sqrt(-(2.*n-2).*log(1-(coh-eps)));
Z    = beta.*(q-beta);
