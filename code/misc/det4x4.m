function [d] = det4x4(x)

%computes determinant of matrix x, using explicit analytic definition if
%size(x,1) == 4, otherwise use matlab det-function

siz = size(x);
if all(siz(1:2)<2)
  error('size of x should be 4x4xsomethingelse');
elseif all(siz(1:2)==4)
  d = x(1,1,:,:).*det3x3(x(2:4,2:4,    :,:)) - ...
      x(1,2,:,:).*det3x3(x(2:4,[1 3 4],:,:)) + ...
      x(1,3,:,:).*det3x3(x(2:4,[1 2 4],:,:)) - ...
      x(1,4,:,:).*det3x3(x(2:4,1:3,    :,:));
% elseif numel(siz)==2
%   d = det(x);
else
  %error   
  %write for loop
  %for
  %end
end
