function wm=winmean(x,percent)

% function wm=winmean(x,percent)
% returns the winsorized mean of x
% x must be a vector
% percent must be between 0 and 100; the function winsorizes the lower and
% upper extreme g values of x, where g=floor((percent/100)*length(x)). If x is empty, then NaN is
% returned. 
% percent=20 by default
%
% Original code provided by Prof. Patrick J. Bennett, McMaster University
% Edit input checks: GAR - University of Glasgow - Nov 2008
% 
% See Rand R. Wilcox (2005), p.27-29
%
% See also WINSAMPLE

if nargin < 2;percent=20;end

% The output size for [] is a special case, handle it here.
if isequal(x,[]), wm = NaN; return; end;

% make sure that x is a vector
sz = size(x);
if sz > 2 | min(sz) > 1    
  error('winsample requires x to be a vector, not a matrix.');
end

wx=winsample(x,percent);
wm=mean(wx);

return
