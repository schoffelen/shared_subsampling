function [sem] = trimsem(x, percent, dim)

%type something here

if nargin<3
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2 || isempty(percent)
  percent = 0.2;
end

n     = size(x, dim);
sigma = winvar(x, percent, dim);
sem   = sqrt(sigma)./((1-2*percent) * sqrt(n));
