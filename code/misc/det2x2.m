function [d] = det2x2(x)

% DET2X2 computes the determinant of a matrix X, using the explicit analytic
% definition if size(X,1) < 4, otherwise use the matlab det()-function. The
% input matrix X can be N-dimensional (up to 4 dimensions), provided the
% first 2 dimensions are equal. This function is also implemented as a
% mex-file. It may be useful when multiple determinants are needed multiple
% times.


siz = size(x);
if all(siz(1:2)==2),
  d = x(1,1,:,:).*x(2,2,:,:) - x(1,2,:,:).*x(2,1,:,:);
elseif all(siz(1:2)==3),
  d = x(1,1,:,:).*x(2,2,:,:).*x(3,3,:,:) - ...
      x(1,1,:,:).*x(2,3,:,:).*x(3,2,:,:) - ...
      x(1,2,:,:).*x(2,1,:,:).*x(3,3,:,:) + ...
      x(1,2,:,:).*x(2,3,:,:).*x(3,1,:,:) + ...
      x(1,3,:,:).*x(2,1,:,:).*x(3,2,:,:) - ...
      x(1,3,:,:).*x(2,2,:,:).*x(3,1,:,:);
elseif numel(siz)==2,
  d = det(x);
else
  %error   
  %write for loop
  %for
  %end
end
