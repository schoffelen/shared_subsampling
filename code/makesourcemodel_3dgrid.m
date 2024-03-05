function leadfield = makesourcemodel_3dgrid(resolution)

if nargin==0
  resolution = 8;
end

load mri
seg = removefields(mri, {'anatomy' 'inside'});
seg.brain = seg.mask; clear mri;
seg = rmfield(seg, 'mask');

cfg = [];
cfg.resolution = resolution;
cfg.unit       = 'mm';
cfg.mri        = seg;
%cfg.smooth     = 'no';
sourcemodel = ft_prepare_sourcemodel(cfg);
  
load grad;
load headmodel;

cfg = [];
cfg.channel = 'MEG';
cfg.sourcemodel = ft_convert_units(sourcemodel, 'm');
cfg.headmodel   = ft_convert_units(headmodel,'m');
cfg.grad        = ft_convert_units(grad,'m');
%cfg.backproject = 'no';
cfg.singleshell.batchsize = 2000;
leadfield = ft_prepare_leadfield(cfg);

insidevec = find(leadfield.inside);
leadfield.v = cell(size(leadfield.leadfield));
for k = 1:numel(insidevec)
  ik  = insidevec(k);
  tmp = leadfield.leadfield{ik};
  [u,s,v] = svd(tmp,'econ');
  leadfield.leadfield{ik}=tmp*v(:,1:2);
  leadfield.v{ik} = v(:,1:2);
end

% inside = reshape(leadfield.inside, leadfield.dim);
% inside(:,:,1:4) = false;
% leadfield.inside = inside(:);

