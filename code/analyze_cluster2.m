function [d, state, num, val, nvox, coeff_v] = analyze_cluster2(clus, sourcemodel, target, sgn)

% ANALYZE_CLUSTER analyzes the array of clusters clus and returns for each element
% the following quantities:
% d     = the distances between the centre of gravity of the legs of the 6D cluster 
%          and the target,
% state = whether the sign of the data within the clusters == sgn
% num   = the number of separate clusters of voxels
% val   = the most extreme value in the cluster
%
% the target pair is determined by the triplet of inputs: dim, inside, and target
% clus can be a single structure, an array of structures, or a cell array of structures
%
% Use as
%  [d,state,num, val] = analyze_cluster(clus, dim, inside, target, sgn)

% Copyright (C) 2008, CCNi Jan-Mathijs Schoffelen
%
% Revamped March 2022
% $Log$

if nargin==3, sgn = 1; end

if iscell(clus)
  nclus  = length(clus);
  d      = cell(1,nclus);
  state  = cell(1,nclus);
  num    = cell(1,nclus);
  val    = cell(1,nclus);
  for k = 1:nclus
    [d{k}, state{k}, num{k}, val{k}, nvox{k}, coeff_v{k}] = analyze_cluster2(clus{k}, sourcemodel, target, sgn);
  end
  return;
end

if length(clus)>1
  nclus  = length(clus);
  %d      = ones(nclus,1)*inf;
  state  = zeros(nclus,1);
  num    = zeros(nclus,1);
  val    = zeros(nclus,1);
  for k = 1:nclus
    [dd, state(k), tmpnum, val(k), nv, cv] = analyze_cluster2(clus(k), sourcemodel, target, sgn);
    num(k, 1:numel(tmpnum)) = tmpnum;
    nvox(k,1:numel(nv))    = nv;
    coeff_v(k,1:numel(cv)) = cv;
    d(k,1:size(dd,2))    = dd;
  end
  return;
end

dim    = sourcemodel.dim;
inside = find(sourcemodel.inside);
pos    = sourcemodel.pos;

if length(clus)==1
  state = sign(clus.val(1))==sgn;
  
  pdim          = prod(dim);
  [indxx,indxy] = ind2sub([pdim pdim], clus.vox);
  
  boolmat = sparse(indxx,indxy,ones(numel(indxx),1),pdim,pdim);
  
  % create a volume in 3D with the 'activated' voxels
  deg   = full(sum(boolmat)./(clus.nvox./2)); % the degree of connectedness of the 'activated' voxels
  cc    = bwconncomp(reshape(deg,dim),26); % cluster in 3D space
  n     = cc.NumObjects;

  nvox  = zeros(1,n);
  coeff_v = zeros(1,n);
  
  d = nan(2,n);

  for m = 1:n
    xm        = cc.PixelIdxList{m};
    nvox(1,m) = numel(xm);
    md        = mean(deg(xm));
    coeff_v(1,m) = sum((deg(xm)-md).^2)./(md.*nvox(1,m));

    pos_ = pos(xm,:);
    if size(pos_, 1)>3
      try
        % fit a convex hull around the points that constitute a local
        % cluster in 3D space
        tri  = convhull(pos_(:,1), pos_(:,2), pos_(:,3));
%         nrm  = normals(pos_-mean(pos_), tri); % check whether the normals are all outward pointing
% 
%         d1 = sum((pos_-mean(pos_)).^2,2);
%         d2 = sum(nrm.^2,2);

        % check whether the target positions are within the hulls
        for mm = 1:size(target,1)
          sa1 = sum(solid_angle(pos_-target(mm,:), tri));
          if abs(sa1)<1000*eps
            % target point 1 is outside the point cloud
            [proj1, dist1] = ptriprojn(pos_(tri(:,1),:), pos_(tri(:,2),:), pos_(tri(:,3),:), target(mm,:), 1);
            d(mm,m) = min(abs(dist1));
          elseif abs(abs(sa1)-4*pi)<1000*eps
            % target point 1 is within the point cloud
            d(mm,m) = 0;
          end
        end
        
      catch 
        % something went wrong with the convex hull creation, default back
        % to the minimum distance
        for mm = 1:size(target,1)
          d(mm,m) = min(sqrt(sum((pos_-target(mm,:)).^2,2)));
        end
      end

    else
      % a convex hull can not be created, so express the distance as the
      % minimum distance of the target points to any of the points in the
      % point cloud
      for mm = 1:size(target,1)
        d(mm,m) = min(sqrt(sum((pos_-target(mm,:)).^2,2)));
      end
    end
  end
  
  for mm = 1:size(target,1)
    [m(mm), idm(mm)] = min(d(mm,:), [], 2);
  end
  
  if numel(m)==2
    % 3 element vector that reflects the minimum distance to the two target
    % dipoles, and with 1/0 whether or not the closets distances are to two
    % different target points
    d = [m(1) m(2) double(isequal([idm(1) idm(2)],[1 2])||isequal([idm(1) idm(2)],[2 1]))];
  else
    d = d(:)';
  end

  if ~isempty(clus.hit)
    num   = 2 - double(sign(clus.hit)==-1);
  else
    num = 1;
  end
  val   = max(abs(clus.val)).*sign(clus.val(1));
else
  d = []; state = []; num = []; val = []; nvox = zeros(0,2); coeff_v = zeros(0,2);
end

function [varargout] = ind2sub(siz,ndx)
%IND2SUB Multiple subscripts from linear index.

siz    = double(siz);
lensiz = length(siz);

k = cumprod(siz);
for i = lensiz:-1:2
  vi = rem(ndx-1, k(i-1)) + 1;
  vj = (ndx - vi)/k(i-1) + 1;
  varargout{i} = double(vj);
  ndx = vi;
end
varargout{1} = double(vi);
