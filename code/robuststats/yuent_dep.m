function [y, mu1, mu2, denom] = yuent_dep(x1, x2, percent, dim, flag)

if nargin<5
  flag = 1;
end

if nargin<4
  dim = find(size(x1)>1, 1, 'first');
end

if nargin<3
  percent = 0.2;
end

n1  =     size(x1, dim);
n2  =     size(x2, dim);

g1  =    floor(percent*n1);
g2  =    floor(percent*n2);

[mu1, thrlo1, thrhi1] = trimmean(x1, percent, dim, flag);
[mu2, thrlo2, thrhi2] = trimmean(x2, percent, dim, flag); 

s1  =   winvar(x1, cat(dim,thrlo1,thrhi1), dim);
s2  =   winvar(x2, cat(dim,thrlo2,thrhi2), dim);

h1  = n1 - 2*g1;
h2  = n2 - 2*g2;

d1  = ( (n1-1)*s1 )/( h1*(h1-1) );
d2  = ( (n2-1)*s2 )/( h2*(h2-1) );

s12 =   wincov(x1, x2, cat(dim,thrlo1,thrhi1), cat(dim,thrlo2,thrhi2), dim);

h12 = h1;
d12 = ( (n1-1)*s12)/( h12*(h12-1) ); 

denom = sqrt(d1+d2-2.*d12);
y     = (mu1-mu2)./denom;
