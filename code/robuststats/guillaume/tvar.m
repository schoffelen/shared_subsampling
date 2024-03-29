function [tv,g]=tvar(x,percent)

% function [tv]=tvar(x,percent)
% Returns the estimated squared standard error of the trimmed mean of x.
% x must be a vector. If x is empty, then NaN is returned.
% percent must be a number between 0 and 100
%
% The trimmed variance is calculated from the winsorized variance, vw,
% according to the formula vw/(k*length(x)), where k=(1-2*percent/100)^2.
%
% Original code provided by Prof. Patrick J. Bennett, McMaster University
%
% See Rand R. Wilcox (2001), Fundamentals of Modern Statisical Methods, page 164.
% See also Rand R. Wilcox (2005), p.61-63
%
% See also WINVAR

% The output size for [] is a special case, handle it here.
if isequal(x,[]), tv = NaN; return; end

if nargin < 2
    error('tvar:TooFewInputs', 'tvar requires two input arguments.');
elseif percent >= 100 || percent < 0
    error('tvar:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

% make sure that x is a vector
sz = size(x);
if sz > 2 | min(sz) > 1    
  error('tvar requires x to be a vector, not a matrix.');
end

[wv,g]=winvar(x,percent);
k=(1-2*percent/100)^2;
tv=wv/(k*length(x));

return
