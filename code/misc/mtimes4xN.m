function [z] = mtimes4xN(x, y)

% MTIMES4XN computes x*y where x = 4x4xM (or Mx4x4) and y = 4xNxM (or
% Mx4xN), and output dimensionatity is 4xNxM (or Mx4xN).
% if M==4 it is treated as Mx4x4

siz   = size(x);
sizy  = [size(y) 1];
z     = complex(zeros(sizy));
%xconj = conj(x);

if siz(2)==siz(3)
  for k = 1:sizy(3)
    z(:,1,k) = x(:,1,1).*y(:,1,k) + x(:,1,2).*y(:,2,k) + x(:,1,3).*y(:,3,k) + x(:,1,4).*y(:,4,k);
    z(:,2,k) = x(:,2,1).*y(:,1,k) + x(:,2,2).*y(:,2,k) + x(:,2,3).*y(:,3,k) + x(:,2,4).*y(:,4,k);
    z(:,3,k) = x(:,3,1).*y(:,1,k) + x(:,3,2).*y(:,2,k) + x(:,3,3).*y(:,3,k) + x(:,3,4).*y(:,4,k);
    z(:,4,k) = x(:,4,1).*y(:,1,k) + x(:,4,2).*y(:,2,k) + x(:,4,3).*y(:,3,k) + x(:,4,4).*y(:,4,k);
  end
elseif siz(1)==siz(2)
  for k = 1:sizy(2)
    z(1,k,:,:) = x(1,1,:,:).*y(1,k,:,:) + x(1,2,:,:).*y(2,k,:,:) + x(1,3,:,:).*y(3,k,:,:) + x(1,4,:,:).*y(4,k,:,:);
    z(2,k,:,:) = x(2,1,:,:).*y(1,k,:,:) + x(2,2,:,:).*y(2,k,:,:) + x(2,3,:,:).*y(3,k,:,:) + x(2,4,:,:).*y(4,k,:,:);
    z(3,k,:,:) = x(3,1,:,:).*y(1,k,:,:) + x(3,2,:,:).*y(2,k,:,:) + x(3,3,:,:).*y(3,k,:,:) + x(3,4,:,:).*y(4,k,:,:);
    z(4,k,:,:) = x(4,1,:,:).*y(1,k,:,:) + x(4,2,:,:).*y(2,k,:,:) + x(4,3,:,:).*y(3,k,:,:) + x(4,4,:,:).*y(4,k,:,:);
  end
else
  error('unsupported dimensionality in the input');
end
