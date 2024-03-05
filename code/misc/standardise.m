function [x,mx,sx] = standardise(x,dim,lim)

% X = STANDARDISE(X, DIM) computes the zscore of a matrix along dimension dim
% has similar functionality as the stats-toolbox's zscore function

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if nargin == 1, 
  dim = find(size(x)>1,1,'first');
  lim = [];
  flag = 0;
elseif nargin == 2,
  lim = [];
  flag = 0;
else
  flag = 1;
end

ndim   = numel(size(x));
if flag
  ix     = cell(1,6);
  for k = 1:numel(ix)
    if k>ndim
      ix{k} = 1;
    else
      ix{k} = 1:size(x,k);
    end
  end
  ix{dim} = lim(1):lim(2);

  mx = mean(x(ix{1},ix{2},ix{3},ix{4},ix{5},ix{6}),dim);
  sx = std(x(ix{1},ix{2},ix{3},ix{4},ix{5},ix{6}),1,dim);
else
  mx = mean(x, dim);
  sx = std(x, 1, dim);
end

switch dim
case 1
  repvec = ones(size(x,1),1);
  x = (x - mx(repvec,:))./sx(repvec,:);
case 2
  repvec = ones(size(x,2),1);
  x = (x - mx(:,repvec))./sx(:,repvec);
otherwise
  repvec = ones(1,ndim);
  repvec(dim) = size(x,dim); 
  x      = (x - repmat(mx,repvec))./repmat(sx,repvec);
end

