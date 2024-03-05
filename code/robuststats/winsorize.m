function [x] = winsorize(x, percent, dim)

% WINSORIZE winsorizes sample x in dimension dim setting the 
% percent percent observations in each of the tails
% to the value of the observation at threshold.
% Implemented according to Wilcox Robust estimation
% and hypothesis testing. If percent is a 2 element vector, it is treated
% as the upper and lower threshold (i.e. a call to percthreshold has
% already been executed).

if nargin<3
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2
  percent = 0.2;
end

switch numel(percent)
  case 1
    [thrlo, thrhi] = percthreshold(x, percent, dim, 0);
  otherwise
    switch dim
      case 1
        thrlo = percent(1,:,:,:);
        thrhi = percent(2,:,:,:);
      case 2
        thrlo = percent(:,1,:,:);
        thrhi = percent(:,2,:,:);
      case 3
        thrlo = percent(:,:,1,:);
        thrhi = percent(:,:,2,:);
      case 4
        thrlo = percent(:,:,:,1);
        thrhi = percent(:,:,:,2);
    end
end
n              = size(x, dim);
%g              = floor(percent*n);
repdimvec      = ones(1, ndims(x));
repdimvec(dim) = n;

%repdimvec(dim) = g;
%x(bsxfun(@lt, x, thrlo)) = repmat(thrlo, repdimvec);
%x(bsxfun(@gt, x, thrhi)) = repmat(thrhi, repdimvec);
x = max(x, repmat(thrlo, repdimvec));
x = min(x, repmat(thrhi, repdimvec));

