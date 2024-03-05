function [z] = sandwichMx4(x, y)

% z = sandwichMx2(x, y)
% compute x*y*x' provided y = hermitian
% and dimensionality is 2x2xN, x should be Mx2xN

%FIXME build in check for hermitianity
z     = complex(zeros(size(x,1),size(x,1),size(x,3)));
xconj = conj(x);
xabs2 = abs(x).^2;

for k = 1:size(x,1)
  z(k,k,:,:) = xabs2(k,1,:,:) .* y(1,1,:,:) + ...
               xabs2(k,2,:,:) .* y(2,2,:,:) + ...
               xabs2(k,3,:,:) .* y(3,3,:,:) + ...
               xabs2(k,4,:,:) .* y(4,4,:,:) + ...
               2.*real( x(k,2,:,:).*y(2,1,:,:).*xconj(k,1,:,:) ) + ...
               2.*real( x(k,3,:,:).*y(3,1,:,:).*xconj(k,1,:,:) ) + ...
               2.*real( x(k,4,:,:).*y(4,1,:,:).*xconj(k,1,:,:) ) + ...
               2.*real( x(k,3,:,:).*y(3,2,:,:).*xconj(k,2,:,:) ) + ...
               2.*real( x(k,4,:,:).*y(4,2,:,:).*xconj(k,2,:,:) ) + ...
               2.*real( x(k,4,:,:).*y(4,3,:,:).*xconj(k,3,:,:) );
    
  for m = (k+1):size(x,1)  
    z(m,k,:,:) = (x(m,1,:,:).*y(1,1,:,:) + x(m,2,:,:).*y(2,1,:,:) + x(m,3,:,:).*y(3,1,:,:) + x(m,4,:,:).*y(4,1,:,:)).*xconj(k,1,:,:) + ...
             (x(m,1,:,:).*y(1,2,:,:) + x(m,2,:,:).*y(2,2,:,:) + x(m,3,:,:).*y(3,2,:,:) + x(m,4,:,:).*y(4,2,:,:)).*xconj(k,2,:,:) + ...
             (x(m,1,:,:).*y(1,3,:,:) + x(m,2,:,:).*y(2,3,:,:) + x(m,3,:,:).*y(3,3,:,:) + x(m,4,:,:).*y(4,3,:,:)).*xconj(k,3,:,:) + ...
             (x(m,1,:,:).*y(1,4,:,:) + x(m,2,:,:).*y(2,4,:,:) + x(m,3,:,:).*y(3,4,:,:) + x(m,4,:,:).*y(4,4,:,:)).*xconj(k,4,:,:);

    
    z(k,m,:,:) = conj(z(m,k,:,:));
  end
end

%b1 b2     a1 a2'   b1' b3'
%b3 b4     a2 a3    b2' b4'
%
%b1*a1+b2*a2  b1*a2'+b2*a3  b1' b3'
%b3*a1+b4*a2  b3*a2'+b4*a3  b2' b4'
%
%b1*a1*b1'+b2*a2*b1'+b1*a2'*b2'+b2*a3*b2' b1*a1*b3'+b2*a2*b3'+b1*a2'*b4'+b2*a3*b4'
%b3*a1*b1'+b4*a2*b1'+b3*a2'*b2'+b4*a3*b2' b3*a1*b3'+b4*a2*b3'+b3*a2'*b4'+b4*a3*b4'
%
%a1*abs(b1)^2 + a2*(b1'*b2) + a2'*(b1*b2') + a3*abs(b2)^2    a1*b1*b3'    + a2*b2*b3'   + a2'*b1*b4'   + a3*b2*b4'
%a1*b1'*b3    + a2*b1'*b4   + a2'*b2'*b3   + a3*b2'*b4       a1*abs(b3)^2 + a2*(b3'*b4) + a2'*(b3*b4') + a3*abs(b4)^2
