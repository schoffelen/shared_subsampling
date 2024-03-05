function [trimmedMean,g]=tm(x,percent)
%
% function [trimmedMean,g]=tm(x,percent)
%
% Returns the trimmed mean of x as trimmedMean.
% returns number of winsorized observations at each tail of sample as g
% x must be a vector
% percent must be a number between 0 and 100
% The trimmed mean is calculated from x after dropping the lower and upper g values from
% x, where g=floor((percent/100)*length(x)). If x is empty, then NaN is
% returned.
%
% See Wilcox (2005), Introduction to Robust Estimation and Hypothesis
% Testing (2nd Edition), Chapter 3, for a description of trimmed means.
%
% Original code provided by Prof. Patrick J. Bennett, McMaster University
%
% See also TSAMPLE

% The output size for [] is a special case, handle it here.
if isequal(x,[]), trimmedMean = NaN; return; end;

% if nargin < 2
%     error('tm:TooFewInputs', 'tm requires two input arguments.');
% elseif percent >= 100 || percent < 0
%     error('tm:InvalidPercent', 'PERCENT must be between 0 and 100.');
% end

% make sure that x is a vector
sz = size(x);
if sz > 2 | min(sz) > 1    
  error('tm requires x to be a vector, not a matrix.');
end

[tx,g]=tsample(x,percent);
trimmedMean=mean(tx);

return
