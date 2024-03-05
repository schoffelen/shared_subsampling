function [v,d] = eig2x2(x)

siz = size(x);
if numel(siz)>3
  error('not yet supported for >3D matrices');
end

if numel(siz)==2
  [v,d] = eig(x);
elseif all(siz(1:2)==2)
  d = zeros(siz);
  v = zeros(siz);
  
  aplusd = x(1,1,:)+x(2,2,:);
  detA   = det2x2(x);
  Discr  = aplusd.^2-4.*detA;
  
  d(1,1,:) = (aplusd+sqrt(Discr))./2;
  d(2,2,:) = (aplusd-sqrt(Discr))./2;

  % if both off-diagonal elements are 0
  b0         = shiftdim(x(1,2,:)==0);
  c0         = shiftdim(x(2,1,:)==0);
  bc0        = b0 & c0;
  v(:,:,bc0) = repmat(eye(2),[1 1 sum(bc0)]);
  
  % if not
%   v(1,1,~bc0) = x(1,2,~bc0);
%   v(1,2,~bc0) = x(1,2,~bc0);
%   v(2,1,~bc0) = d(1,~bc0) - shiftdim(x(1,1,~bc0),1);
%   v(2,2,~bc0) = d(2,~bc0) - shiftdim(x(1,1,~bc0),1);
  
  v(1,1,~c0) = d(1,1,~c0) - x(2,2,~c0);
  v(1,2,~c0) = d(2,2,~c0) - x(2,2,~c0);
  v(2,1,~c0) = x(2,1,~c0);%d(1,~bc0);% - shiftdim(x(1,1,~bc0),1);
  v(2,2,~c0) = x(2,1,~c0);%d(2,~bc0);% - shiftdim(x(1,1,~bc0),1);
  
  % norm normalize
  v = v./repmat(sqrt(sum(v.^2,1)), [2 1 1]);
  
end