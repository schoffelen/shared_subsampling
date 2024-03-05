function [mu, thrlo, thrhi] = trimmean(x, percent, dim, flag)

%Replacement for the matlab version of trimmean
%Implementation according to Wilcox; Robust Estimation
%and Hypothesis Testing
%
%The difference with the matlab implementation is that
%percent is applied to both tails of the distribution
%in the present implementation, whereas the matlab
%implementation divides percent by 2 first.
%
%Use as [MU] = TRIMMEAN(X, PERCENT, DIM)
%
%This function gives output which is different from matlab

if nargin<4
  flag = 1;
end

if nargin<3
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2 || isempty(percent)
  percent = 0.2;
end

n              = size(x, dim);
g              = floor(percent*n);

if flag
  [thrlo, thrhi, ix] = percthreshold(x, percent, dim, 2);
  
  % the following deals correctly with ties, but is much slower.
  % when working with real data the probability to get a tie is low I would
  % think.
  denom = n-2*g;
  mu = sum(x.*double(ix),dim)./denom;
  %siz = size(x);
  %siz(dim) = 1;
  %mu  = zeros(siz);
  % switch dim
  %   case 1
  %     for k=1:size(x,2)
  %       mu(:,k,:,:,:,:,:) = sum(x(ix(:,k),k,:,:,:,:))./denom;
  %     end
  %   case 2
  %     mu = sum(x.*double(ix),2)./denom;
  %     %for k=1:size(x,1)
  %     %  mu(k,:,:,:,:,:) = sum(x(k,ix(k,:),:,:,:,:,:),2)./denom;
  %     %end
  %   case 3
  %     for k=1:size(x,1)
  %       for m=1:size(x,2)
  %         mu(k,m,:,:,:,:) = sum(x(k,m,ix(k,m,:),:,:,:),3)./denom;
  %       end
  %     end
  %   case 4
  %     for k=1:size(x,1)
  %       for m=1:size(x,2)
  %         for n=1:size(x,3)
  %           mu(k,m,n,:,:,:,:) = sum(x(k,m,n,ix(k,m,n,:),:,:,:),4)./denom;
  %         end
  %       end
  %     end
  %   otherwise
  %     error('trimmean for dim>4 not yet supported');
  % end
else
  
  [thrlo, thrhi] = percthreshold(x, percent, dim, 0);
  
  % the bsxfun does not deal correctly with ties, not sure whether this makes
  % a difference in practice
  x(bsxfun(@lt, x, thrlo) | bsxfun(@gt, x, thrhi)) = 0;
  mu = sum(x, dim)./(n-2*g);
end