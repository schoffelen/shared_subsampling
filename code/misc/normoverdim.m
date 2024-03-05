function [y, x_norm] = normoverdim(x, dim)

x_norm = sqrt(sum(abs(x).^2,dim));
y      = bsxfun(@rdivide, x, x_norm);
