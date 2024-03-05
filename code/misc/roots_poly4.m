function out = roots_poly4(c,flag)

if nargin<2
  flag = false;
end


delta0 = c(:,3).^2 - 3.*c(:,2).*c(:,4) + 12.*c(:,5);
delta1 = 2.*(c(:,3).^3) - 9.*c(:,2).*c(:,3).*c(:,4) + ...
         27.*(c(:,2).^2).*c(:,5) + 27.*c(:,4).^2 - ...
         72.*c(:,3).*c(:,5);

p = c(:,3) - (3./8).*(c(:,2).^2);
q = (c(:,2).^3)./8 -(c(:,2).*c(:,3))./2 + c(:,4);

Q = ((delta1 + sqrt(delta1.^2 - 4.*(delta0.^3)))./2).^(1/3);
S = sqrt(-(2./3).*p + (Q+delta0./Q)./3)./2;

tmp0 = -c(:,2)./4 - S;
tmp1 = sqrt( - 4.*S.^2 - 2.*p + q./S )./2;

if ~flag
  % c = Nx5, with leading coefficient assumed to be 1
  out = zeros(size(c,1),4);
  
  out(:,1) = tmp0 - tmp1;
  out(:,2) = tmp0 + tmp1;
  
  tmp0 = -c(:,2)./4 + S;
  tmp1 = sqrt( - 4.*S.^2 - 2.*p - q./S )./2;
  
  out(:,3) = tmp0 - tmp1;
  out(:,4) = tmp0 + tmp1;
  
  out = real(out);
else
  out = zeros(size(c,1),1);
  
  tmp0 = -c(:,2)./4 + S;
  tmp1 = sqrt( - 4.*S.^2 - 2.*p - q./S )./2;
  
  out(:) = real(tmp0 + tmp1);
  
end
