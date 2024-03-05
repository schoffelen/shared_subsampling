function [a, r] = fitcubic(x, y)

% fit a cubic with only a non-zero 3d and 1st power (because its symmetry
% axis should be the y-axis, and it should go through [0, 0]

X = [x(:) x(:).^3];
a = pinv(X)*y;

yhat = X*a;

r = sum((y-yhat).^2)./sum(y.^2);
