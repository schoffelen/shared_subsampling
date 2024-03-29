function [wx,g]=winsample(x,percent)

% function [wx,g]=winsample(x,percent)
% returns the winsorized sample of x as wx
% returns number of winsorized observations at each tail of sample as g
% x must be a vector
% percent must be between 0 and 100; the function converts the lower and
% upper extreme g values of x to x(g+1) and x(n-g), respectively, where
% g=floor((percent/100)*length(x)). If x is empty, then NaN is returned. 
%
% Original code provided by Prof. Patrick J. Bennett, McMaster University
%
% Modified to avoid sorting the output: GAR - University of Glasgow - Dec 2007
% Edit input checks: GAR - University of Glasgow - Nov 2008
%
% See also WINVAR

if nargin < 2;percent=20;end

% The output size for [] is a special case, handle it here.
if isequal(x,[]), wx = NaN; return; end;

% make sure that x is a vector
sz = size(x);
if sz > 2 | min(sz) > 1    
  error('winsample requires x to be a vector, not a matrix.');
end

n = length(x);
xsort=sort(x);
g=floor((percent/100)*n);

loval=xsort(g+1);
hival=xsort(n-g);

% wx=xsort;
% 
% wx(1:g)=loval+zeros(1,g);
% wx(n-g+1:n)=hival+zeros(1,g);

wx = x;
wx(find(wx<=loval)) = loval;
wx(find(wx>=hival)) = hival;

return
