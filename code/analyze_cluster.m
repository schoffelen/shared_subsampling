function [d, state, num, val, mindist] = analyze_cluster(clus, dim, inside, targetindx, sgn)

% ANALYZE_CLUSTER analyzes the array of clusters clus and returns for each element
% the following quantities:
% d     = the euclidean distance between the pair of voxels closest to the target pair,
%        and the target pair
% state = whether the sign of the data within the clusters == sgn
% num   = the number of separate clusters of voxels
% val   = the most extreme value in the cluster
%
% the target pair is determined by the triplet of inputs: dim, inside, and targetindx
% clus can be a single structure, an array of structures, or a cell array of structures
%
% Use as
%  [d,state,num, val] = analyze_cluster(clus, dim, inside, targetindx, sgn)

% Copyright (C) 2008, CCNi Jan-Mathijs Schoffelen
% $Log$

if nargin==4, sgn = 1; end

if iscell(clus)
  nclus  = length(clus);
  d      = cell(1,nclus);
  state  = cell(1,nclus);
  num    = cell(1,nclus);
  val    = cell(1,nclus);
  for k = 1:nclus
    [d{k}, state{k}, num{k}, val{k}, mindist{k}] = analyze_cluster(clus{k}, dim, inside, targetindx, sgn);
  end
  return
end

if length(clus)>1
  nclus  = length(clus);
  d      = ones(nclus,1)*inf;
  state  = zeros(nclus,1);
  num    = zeros(nclus,1);
  val    = zeros(nclus,1);
  mindist = zeros(nclus,1);
  for k = 1:nclus
    [d(k), state(k), tmpnum, val(k), mindist(k)] = analyze_cluster(clus(k), dim, inside, targetindx, sgn);
    num(k, 1:numel(tmpnum)) = tmpnum;
  end
  return
end

if length(clus)==1
  ix        = [];
  [ix(:,1),ix(:,2),ix(:,3),ix(:,4),ix(:,5),ix(:,6)] = ...
       ind2sub([dim dim], clus.vox); % vox is actually a misnomer, because it refers to edges.
  
  d = inf;
  for kk = 1:size(targetindx,3)
    dpos  = ix - repmat([targetindx(1,:,kk) targetindx(2,:,kk)], [size(ix,1) 1]);
    mind  = sqrt(sum(dpos.^2,2)); % this is apparently a 6D distance, expressed in units of grid spacing
    dpos  = ix - repmat([targetindx(2,:,kk) targetindx(1,:,kk)], [size(ix,1) 1]);
    mind  = min(min(mind,sqrt(sum(dpos.^2,2))));
    d     = min(mind, d);
  end
  state = sign(clus.val(1))==sgn;
  num   = 2 - double(sign(clus.hit)==-1);
  val   = max(abs(clus.val)).*sign(clus.val(1));
  
  vox1 = unique(ix(:,1:3),'rows'); n1 = size(vox1,1);
  dum  = zeros(dim);
  for k = 1:n1
    dum(vox1(k,1),vox1(k,2),vox1(k,3)) = 1;
  end
  [dum, nclus] = bwlabeln(dum, 26);
  if nclus>1
    mindist = 0;
  else
    mindist = 1;
  end  
else
  d = []; state = []; num = []; val = []; mindist = [];
end
