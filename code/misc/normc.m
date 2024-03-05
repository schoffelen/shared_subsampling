function [y, x_norm] = normc(x)

% NORMC column-wise norm normalisation

x_norm = sqrt(sum(x.^2));
y      = full(x*spdiags(1./x_norm(:),0,size(x,2),size(x,2)));
