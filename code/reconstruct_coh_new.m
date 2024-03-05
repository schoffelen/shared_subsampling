function [coh, C, Cc, Corig, Corig1, c, c0, hl] = reconstruct_coh_new(freq, freqnoise, leadfield, opts)

if nargin<4
  opts = [];
end

N       = ft_getopt(opts, 'N',      100); % number of sensors in the subsample
nrand   = ft_getopt(opts, 'nrand',  0);
refindx = ft_getopt(opts, 'refindx', []);
lambda  = ft_getopt(opts, 'lambda', 0);
lambda_sub = ft_getopt(opts, 'lambda_sub', 0);
fixedori   = ft_getopt(opts, 'fixedori', 'lambda1');
outputpow  = ft_getopt(opts, 'outputpow', false);
state      = ft_getopt(opts, 'state', 'shuffle');
batchsize  = ft_getopt(opts, 'batchsize', 100);
invTol     = ft_getopt(opts, 'invTol', 1e-12);
fitmethod  = ft_getopt(opts, 'fitmethod', 1);
avgflag    = ft_getopt(opts, 'avgflag', 1);
outputflags = ft_getopt(opts, 'outputflags', [1 1 0 0 0]); % outputflags for the subsampling
ndip_int    = ft_getopt(opts, 'ndip_int', 2); % old default is 2
subsampletrials = ft_getopt(opts, 'subsampletrials', 1);
lfproject  = ft_getopt(opts, 'lfproject2two_before', 1); % flag that determines whether the leadfield will be projected onto the plane prior to subsample, or after (default is prior).
% note: the lfproject==false only makes sense if leadfields have been computed with rank3

subsample.N     = N;
subsample.nrand = nrand;


[a,b]                                 = match_str(leadfield.label, freq.label);
leadfield.leadfield(leadfield.inside) = cellrowselect(leadfield.leadfield(leadfield.inside),a);
leadfield                             = removefields(leadfield, 'v');

if ~isempty(refindx)
  coh     = estimate_coh2x2_2dip_subsample(leadfield,freq,freqnoise,'lambda',lambda,    'memory','high','outputflags',[1 1 1 1 1],'refindx',refindx, 'fixedori', fixedori, 'outputpow', outputpow, 'invTol', invTol, 'subsampletrials', subsampletrials);
  coh_sub = estimate_coh2x2_2dip_subsample(leadfield,freq,freqnoise,'lambda',lambda_sub,'memory','high','outputflags',outputflags,'refindx',refindx, 'fixedori', fixedori, 'outputpow', outputpow, 'subsample', subsample, 'state', state, 'invTol', invTol, 'fitmethod', fitmethod, 'ori0', coh.ori, 'avgflag', avgflag, 'subsampletrials', subsampletrials, 'lfproject2two_before', lfproject);
else
  nbatch    = ceil(sum(leadfield.inside)/batchsize);
  for k = 1:nbatch
    this_refindx = (k-1)*batchsize + (1:batchsize);
    this_refindx(this_refindx>sum(leadfield.inside)) = [];
    coh(k)     = estimate_coh2x2_2dip_subsample(leadfield,freq,freqnoise,'lambda',lambda,    'memory','high','outputflags',[1 1 1 1 1],'refindx',this_refindx, 'fixedori', fixedori, 'outputpow', outputpow, 'invTol', invTol, 'subsampletrials', subsampletrials);
    coh_sub(k) = estimate_coh2x2_2dip_subsample(leadfield,freq,freqnoise,'lambda',lambda_sub,'memory','high','outputflags',outputflags,'refindx',this_refindx, 'fixedori', fixedori, 'outputpow', outputpow, 'subsample', subsample, 'state', state, 'invTol', invTol, 'fitmethod', fitmethod, 'ori0', coh(k).ori, 'avgflag', avgflag, 'subsampletrials', subsampletrials, 'lfproject2two_before', lfproject);
  end
end

% sanity check
ref = cat(2,coh(:).refindx);
if ~isequal(ref,1:sum(leadfield.inside))
  C = coh_sub;
  Cc = [];
  Corig = [];
  Corig1 = [];
  try
    % try to output the subsampled 2-dipole results
    c  = abs(cat(2,coh_sub(:).coh)./cat(2,coh_sub(:).dvar))./sqrt(numel(coh_sub(1).nsub));
    c0 = abs(cat(2,coh_sub(:).coh0)./cat(2,coh_sub(:).dvar))./sqrt(numel(coh_sub(1).nsub));
  catch
    try
      % if the above fails, try to output the 1-dipole results
      c  = abs(cat(2,coh_sub(:).coh1)./cat(2,coh_sub(:).dvar1))./sqrt(numel(coh_sub(1).nsub));
      c0 = abs(cat(2,coh_sub(:).coh1_0)./cat(2,coh_sub(:).dvar1))./sqrt(numel(coh_sub(1).nsub));
    catch
      c  = [];
      c0 = [];
    end
  end
  if isfield(opts, 'ix') && ~isempty(opts.ix)
    [hl(:,1),hl(:,2),hl(:,3)]=ind2sub(leadfield.dim,opts.ix);
  end
  
else
  
  if all(isfinite(opts.ix))
    cfg        = [];
    cfg.dim    = leadfield.dim;
    cfg.inside = find(leadfield.inside);

    % make an Nx2 indx matrix that specifies the indices of the
    % interactions, one edge per row
    tmp1 = repmat(1:ndip_int, ndip_int, 1);
    tmp2 = repmat((1:ndip_int)', 1, ndip_int);
    tmp1 = tmp1(tril(ones(ndip_int),-1)>0);
    tmp2 = tmp2(tril(ones(ndip_int),-1)>0);

    cfg.indx   = opts.ix([tmp1(:) tmp2(:)]);
    cfg.tail   = 1;

    dcoh = abs(cat(2,coh(:).coh))-abs(cat(2,coh(:).coh0));
    dcoh = (dcoh+dcoh')./2;
    thr  = prctile(dcoh(:),[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]);
    Corig = cell(1,numel(thr));
    for k = 1:numel(thr)
      cfg.threshold = thr(k);
      Corig{k} = cluster6D(cfg,dcoh);
    end

    dcoh = abs(cat(2,coh(:).coh1))-abs(cat(2,coh(:).coh1_0));
    dcoh = (dcoh+dcoh')./2;
    thr  = prctile(dcoh(:),[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]);
    Corig1 = cell(1,numel(thr));
    for k = 1:numel(thr)
      cfg.threshold = thr(k);
      Corig1{k} = cluster6D(cfg,dcoh);
    end

    dcoh = abs(cat(2,coh(:).cohc));
    dcoh = (dcoh+dcoh')./2;
    thr  = prctile(dcoh(:),[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]);
    Cc   = cell(1,numel(thr));
    for k = 1:numel(thr)
      cfg.threshold = thr(k);
      Cc{k} = cluster6D(cfg,dcoh);
    end

    dcoh = ((abs(cat(2,coh_sub(:).coh))-abs(cat(2,coh_sub(:).coh0)))./cat(2,coh_sub(:).dvar))./sqrt(numel(coh_sub(1).nsub));
    dcoh = (dcoh+dcoh')./2;
    thr  = prctile(dcoh(:),[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]);
    C = cell(1,numel(thr));
    for k = 1:numel(thr)
      cfg.threshold = thr(k);
      C{k} = cluster6D(cfg,dcoh);
    end

    if ~isempty(opts.ix)
      [hl(:,1),hl(:,2),hl(:,3)]=ind2sub(leadfield.dim,opts.ix);
    end
  end

  tmp = coh;
  clear coh;
  coh = tmp(1);
  coh.coh    = cat(2,tmp(:).coh);
  coh.coh0   = cat(2,tmp(:).coh0);
  coh.cohc   = cat(2,tmp(:).cohc);
  coh.coh1   = cat(2,tmp(:).coh1);
  coh.coh1_0 = cat(2,tmp(:).coh1_0);
  coh.cohsub = cat(2,coh_sub(:).coh);
  coh.cohsub0 = cat(2,coh_sub(:).coh0);
  try, coh.dvar   = cat(2,coh_sub(:).dvar); end
  coh.a      = cat(1,tmp(:).a);
  coh.r      = cat(1,tmp(:).r);
  coh.a1     = cat(1,tmp(:).a1);
  coh.r1     = cat(1,tmp(:).r1);
  coh.refindx = cat(2,tmp(:).refindx);
  clear tmp;
  
  c = abs(cat(2,coh_sub(:).coh)./cat(2,coh_sub(:).dvar))./sqrt(numel(coh_sub(1).nsub));
  c0 = abs(cat(2,coh_sub(:).coh0)./cat(2,coh_sub(:).dvar))./sqrt(numel(coh_sub(1).nsub));
end
