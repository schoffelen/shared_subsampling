function [clus] = cluster6D(cfg, dat)

% [CLUS] = CLUSTER6D(CFG, DAT)
%
% Performs a 6-dimensional clustering on a thresholded data matrix.
% The data matrix is a 2-dimensional representation of a 6-dimensional volume
% Required fields for the input configuration are:
%  cfg.inside
%  cfg.dim
%  cfg.threshold

% Copyright (C) 2008, Jan-Mathijs Schoffelen
% $Log: cluster6D.m,v $
% Revision 1.2  2008/09/09 09:33:43  jan
% extensive update of the original function, leading to new functionality
%

%diagonal is correctly dealt with when two-tailed test. if tail~=0 diagonal has te be manually adjusted

cfg = ft_checkconfig(cfg, 'required', {'inside', 'dim', 'threshold'});

cfg.indx      = ft_getopt(cfg, 'indx',      []);
cfg.conn      = ft_getopt(cfg, 'conn',      26);
cfg.tail      = ft_getopt(cfg, 'tail',      0);
cfg.targetsgn = ft_getopt(cfg, 'targetsgn', 1);
cfg.targetpos = ft_getopt(cfg, 'targetpos', []);
cfg.method    = ft_getopt(cfg, 'method',    'fixed');
cfg.maskdist  = ft_getopt(cfg, 'maskdist',  0);
cfg.threshold2 = ft_getopt(cfg, 'threshold2', 0);
cfg.keepsinglecluster = ft_getopt(cfg, 'keepsinglecluster', false); % also return the suprathreshold blobs which are spatially contiguous

hasindx   = ~isempty(cfg.indx);
hastrgpos = ~isempty(cfg.targetpos);
if hastrgpos
  if ~isfield(cfg, 'pos'), error('cfg.targetpos should be accompanied by cfg.pos'); end
end


%% call recursively if two-tailed
%--------------------------------
if cfg.tail==0
  nvox = size(dat,1);
  if length(cfg.threshold)==1
    t                = abs(cfg.threshold);
    fprintf('running a two-tailed clustering with thresholds %d and %d\n', -t, t);
    tmpcfg           = cfg;
    if strcmp(cfg.method, 'fixed')
      tmpcfg.tail      = 1;
      tmpcfg.threshold = t;
      tmpcfg.indx      = cfg.indx;
      tmpdat           = dat - diag(diag(dat)) + eye(nvox).*tmpcfg.threshold.*1.5;
      clus1            = cluster6D(tmpcfg, tmpdat);
      tmpcfg.tail      = -1;
      tmpcfg.threshold = -t;
      clus2            = cluster6D(tmpcfg, tmpdat);
    elseif strcmp(cfg.method, 'dynamic')
      tmpdat = dat;
      
      tmpcfg.tail      = 1;
      tmpcfg.threshold = t;
      tmpcfg.threshold2 = cfg.threshold2;
      clus1            = cluster6D(tmpcfg, tmpdat);
      tmpcfg.tail      = -1;
      tmpcfg.threshold = -t;
      tmpcfg.threshold2 = -cfg.threshold2;
      clus2            = cluster6D(tmpcfg, tmpdat);
    end
    if ~isempty(clus1) && ~isempty(clus2),
      clus = cat(2, clus1, clus2);
    elseif isempty(clus1),
      clus = clus2;
    elseif isempty(clus2),
      clus = clus1;
    end
  elseif length(cfg.threshold)>1
    t                = cfg.threshold;
    fprintf('running a two-tailed clustering with thresholds %d and %d\n', t(1), t(2));
    %interpret threshold as a lower and upper critical value
    tmpcfg           = cfg;
    tmpcfg.tail      = 1;
    tmpcfg.threshold = cfg.threshold(2);
    clus1            = cluster6D(tmpcfg, tmpdat);
    tmpcfg.tail      = -1;
    tmpcfg.threshold = cfg.threshold(1);
    clus2            = cluster6D(tmpcfg, tmpdat);
    clus             = cat(2, clus1, clus2);
  end
  return;
end
%-------------------------------------------

%% the real computations start here

% initialize some variables
nvox      = size(dat,1);
dim       = cfg.dim;
inside    = cfg.inside(:)';
threshold = cfg.threshold;
if hasindx
  indx = cfg.indx + (cfg.indx(:, [2 1])-1).*prod(dim);
end

% ensure that the diagonal is exceeding the threshold
dat     = dat - diag(diag(dat)) + eye(nvox).*threshold.*1.5;

% perform columnwise clustering of thresholded 3D-data
tmp     = zeros(dim);
clusmat = zeros(prod(dim),'int16');
cnt     = 0;
msk     = -cfg.maskdist:cfg.maskdist;
for k = 1:length(inside)
  if mod(k,100)==0, fprintf('thresholding the data for seed voxel %d/%d\n',k,length(inside)); end
  % remove voxels connected to diagonal
  switch lower(cfg.method)
    case 'fixed'
      if cfg.tail==1
        tmp(inside) = dat(:,k)>threshold;
      elseif cfg.tail==-1
        tmp(inside) = dat(:,k)<threshold;
      end
    case 'dynamic'
      tmp(inside) = dat(:,k); tmp(~isfinite(tmp)) = 0;
      tmp         = 2.*tmp-imhmax(tmp, cfg.threshold)-imhmin(tmp, cfg.threshold);
      if cfg.tail==1
        tmp(inside)  = (dat(:,k) > cfg.threshold2);
      elseif cfg.tail==-1
        tmp(inside)  = (dat(:,k) < cfg.threshold2);
      end
  end
  [ix,iy,iz]     = ind2sub(size(tmp),inside(k));
  ix = ix+msk;ix = min(ix,size(tmp,1));ix = max(1,ix);
  iy = iy+msk;iy = min(iy,size(tmp,2));iy = max(1,iy);
  iz = iz+msk;iz = min(iz,size(tmp,3));iz = max(1,iz);
  tmp(ix,iy,iz)  = true;
  
  tmp2 = bwlabeln(tmp, cfg.conn);
  % put clusters connected to the 'diagonal' to zero.
  % this makes the matrix asymmetric!
  tmp2(tmp2==tmp2(inside(k))) = 0;
  
  % set back to false: this is necessary for the next iteration!
  tmp(ix,iy,iz) = false;
  clusnum       = 1:max(tmp2(:));
  for m = 1:length(clusnum)
    sel = find(tmp2==clusnum(m));
    if ~isempty(sel)
      cnt = cnt + 1;
      clusmat(sel, inside(k)) = cnt;
    end
  end
end
clear ix iy iz;

% threshold
clusmat = clusmat>0;

% make symmetric again and convert to 6-dimensional
clusmat = reshape(clusmat & clusmat', [dim dim]);

% make clusters in 6-dimensional space
clusmat = reshape(bwlabeln(clusmat, conndef(6, 'max')), [prod(dim) prod(dim)]); clear tmp
%clusmat = reshape(bwlabeln_blocked(tmp, conndef(6, 'min'), 6), [prod(dim) prod(dim)]); clear tmp

clusmats = sparse(clusmat);
nclus   = full(max(clusmats(:)));
clusind = 1:nclus;
cnt     = 0;
clus    = struct('vox',[],'val',[],'hit',[],'nvox',[], 'thr', []);

% create the output structure after post-processing of the clusmat-matrix
if 1%nclus<1500,
  while ~isempty(clusind)
    fprintf('processing cluster, %d ', clusind(1));
    [x, y] = find(clusmats==clusind(1));
    nvox   = length(x);
    if nvox>0
      indvec2 = y+prod(dim)*(x-1);
      val     = full(clusmats(indvec2));
      if mode(val)==clusind(1)
        if all(val==val(1)) && sum(val==val(1))==nvox
          fprintf('consisting of one single connected group of voxels\n');
          % cluster consists of one connected group of voxels
          cnt              = cnt+1;
          indvec1          = x+prod(dim)*(y-1);
          clus(cnt).vox    = int32(unique([indvec1;indvec2]));
          clus(cnt).nvox   = length(clus(cnt).vox);
          if hasindx 
            for kk = 1:size(indx,1)
            clus(cnt).hit(1,kk) = -1 - double(sum(ismember(clus(cnt).vox, indx(kk,:)'))==2); 
            end
          end
          clus(cnt).thr    = threshold;
        end
      elseif all(val==val(1)) && sum(val==val(1))==nvox
        fprintf('which is the same cluster as cluster %d\n', val(1));
        % or two clusters belong together with x and y swapped
        cnt              = cnt+1;
        indvec1          = x+prod(dim)*(y-1);
        clus(cnt).vox    = int32(unique([indvec1;indvec2]));
        clus(cnt).nvox   = length(clus(cnt).vox);
        if hasindx
          for kk = 1:size(indx,1)
            clus(cnt).hit(1,kk) = double(sum(ismember(clus(cnt).vox, indx(kk,:)'))==2); 
          end
        end
        clus(cnt).thr    = threshold;
      else
        % don't know how to interpret this
      end
    end
    clusind = sort(setdiff(clusind, [clusind(1) val(1)]));
  end
  
  nvox    = [clus.nvox];
  [~,ind] = sort(nvox, 'descend');
  clus    = clus(ind);
  
  clusmat(inside, inside) = dat;
  for k = 1:length(clus)
    clus(k).val = single(clusmat(clus(k).vox));
  end
else
  %fprintf('there are  %d distinct clusters\n', nclus);
  clus = repmat(clus, [1 nclus]);
end

if ~isempty(cfg.indx)
  for kk = 1:size(indx,1)
    [ix(:,1,kk), ix(:,2,kk), ix(:,3,kk)] = ind2sub(dim, cfg.indx(kk,:));
  end
  [d, state, num, ~, mindist] = analyze_cluster(clus, dim, inside, ix, cfg.targetsgn);
  for k = 1:length(clus)
    clus(k).d     = d(k);
    clus(k).state = state(k);
    clus(k).num   = num(k,:);
    clus(k).mindist = mindist(k);
  end
  if ~isempty(clus) && ~cfg.keepsinglecluster
    clus = clus([clus.mindist]==0);
  end
end
