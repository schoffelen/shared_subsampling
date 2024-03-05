function [y] = yuent(x1, x2, percent, dim)

if nargin<4
  dim = find(size(x1)>1, 1, 'first');
end

if nargin<3
  percent = 0.2;
end

n1  =     size(x1, dim);
g1  =    floor(percent*n1);
[mu1,thrlo1,thrhi1] = trimmean(x1, percent, dim);
s1  =   winvar(x1, cat(dim, thrlo1, thrhi1), dim);

h1  = n1 - 2*g1;
d1  = ( (n1-1)*s1 )/( h1*(h1-1) );

n2  =     size(x2, dim);
g2  =    floor(percent*n2);
[mu2,thrlo2,thrhi2] = trimmean(x1, percent, dim);
s2  =   winvar(x2, cat(dim, thrlo2, thrhi2), dim);

h2  = n2 - 2*g2;
d2  = ( (n2-1)*s2 )/( h2*(h2-1) );

y   = (mu1-mu2)./sqrt(d1+d2);
