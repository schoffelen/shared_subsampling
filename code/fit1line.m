%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use linear regression to fit the 'slope', i.e. noise scaling needed to
% get the coh0 and coh aligned
function [a,r,m,n] = fitslope(x,y)

xorig = x;

x = abs(x);
y = abs(y);

%[t1, t2] = percthreshold(x, 0.5, 1, 0);
%t1  = prctile(x, 95);

if 1
  t1  = median(x);%prctile(x, 50);
  sel = x>t1.*0.5;
  
  x  = x(sel);
  mx = mean(x);
  
  y  = y(sel);
  my = mean(y);
  
  %w = spdiags(x(:), 0, n, n);
  x = x-mx;
  y = y-my;
  %a = (x'*w*y)./(x'*w*x);
  xsq = x.^2;
  a = (xsq'*y)./(xsq'*x); % this is the same as the weighted regression above, with weights the values of x
  r = sum((y-a*x).^2)./(y'*y);
  
  if nargout>2
    m = a*xorig;
    n = sum(sel);
  end
else
  [srt, ix] = sort(x);
  
  krn = ones(999,1)./999;
  xm_(ix,1) = convn(x(ix(:)),krn,'same');
  ym_(ix,1) = convn(y(ix(:)),krn,'same');
  
  x = x-xm_;
  y = y-ym_;
  
  xdenom(ix,1) = convn(x(ix(:)).^2,krn,'same');
  xy(ix,1)     = convn(x(ix(:)).*y(ix(:)),krn,'same');
  a(:,1)       = xy./xdenom; 
  r            = sum((y-a.*x).^2)./(y'*y);
  
end
