function [wv,g]=winvar(x,percent)

% function [wv,g]=winvar(x,percent)
% returns the winsorized variance of x as wv
% returns number of winsorized observations at each tail of sample as g
% x must be a vector
% percent must be between 0 and 100; the function winsorizes the lower and
% upper extreme g values of x, where g=floor((percent/100)*length(x)).  If x is empty, then NaN is
% returned. 
%
% Original code provided by Prof. Patrick J. Bennett, McMaster University
% Edit input checks: GAR - University of Glasgow - Nov 2008
%
% See also WINSAMPLE

if nargin < 2;percent=20;end

% The output size for [] is a special case, handle it here.
if isequal(x,[]), wv = NaN; return; end

if nargin < 2
    error('winvar:TooFewInputs', 'winvar requires two input arguments.');
elseif percent >= 100 || percent < 0
    error('winvar:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

% make sure that x is a vector
sz = size(x);
if sz > 2 | min(sz) > 1    
  error('winvar requires x to be a vector, not a matrix.');
end

[wx,g]=winsample(x,percent);
wv=var(wx);

return
