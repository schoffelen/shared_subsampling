function [freq, freqnoise, params] = simulate_data(opts)

ndip    = ft_getopt(opts, 'ndip',   50);
phs     = ft_getopt(opts, 'phs',    pi/4);
amp     = ft_getopt(opts, 'amp',    0.4); % magnitude of correlation
snr     = ft_getopt(opts, 'snr',    0.5); % SNR parameter for simulation sensor noise versus emptyroom signal
snr2    = ft_getopt(opts, 'snr2',   0.9); % SNR parameter for amplitude of sources of interest versus brain noise sources
rho     = ft_getopt(opts, 'rho', []);
ampl    = ft_getopt(opts, 'ampl', []);
patchsigma = ft_getopt(opts, 'patchsigma', 0.4);
noise   = ft_getopt(opts, 'noise', 'emptyroom');

load(sprintf('candidates_patch_%03d',round(100*patchsigma)));
load leadfield;
leadfield = removefields(leadfield, 'v');

if strcmp(noise, 'emptyroom')
  load noisecov
  % remove 'trials' 84 and 88 from freq, they are outliers in terms of the
  % csd
  tmpcfg = [];
  tmpcfg.trials = setdiff(1:numel(freq.cumtapcnt), [84 88]);
  freq = ft_selectdata(tmpcfg, freq);
elseif strcmp(noise, 'task')
  load noisecov_taskdata

  %tmpcfg = [];
  %tmpcfg.trials = 101:200; % always take the same trials
  %freq = ft_selectdata(tmpcfg, freq);
else
  ft_error('the noise option should be either emptyroom or task')
end

% split the 100 trials into 2 times 50, to have non-overlapping data for
% the 'noise condition', and the additive noise for the 'active condition'
rptsel = randperm(numel(freq.cumtapcnt));
rptsel1 = sort(rptsel(1:floor(numel(rptsel)/2)));
rptsel2 = sort(rptsel(ceil(numel(rptsel)/2):end));

[a,b] = match_str(freq.label,label);
tmpcfg = [];
tmpcfg.channel = freq.label(a);
tmpcfg.frequency = 10;
tmpcfg.trials    = rptsel1;
freqnoise = ft_selectdata(tmpcfg,freq);
tmpcfg.trials    = rptsel2;
freqnoise2 = ft_selectdata(tmpcfg, freq); % this is used to simulate the 'active' data

%freqnoise = ft_checkdata(freqnoise,'cmbstyle','fullfast');
%freqnoise.crsspctrm = freqnoise.crsspctrm./(trace(freqnoise.crsspctrm)./size(freqnoise.crsspctrm,1));
%freqnoise.crsspctrm = freqnoise.crsspctrm./norm(freqnoise.crsspctrm,'fro');
%freqnoise.crsspctrm = (1-snr).*freqnoise.crsspctrm;

tmp   = ft_checkdata(freqnoise, 'cmbstyle', 'fullfast');
denom = sqrt(norm(tmp.crsspctrm, 'fro'));
freqnoise.fourierspctrm = (1-snr).*freqnoise.fourierspctrm./denom; % scale with noise parameter, to stay consistent with the 'active data'

%freqnoise2 = ft_checkdata(freqnoise2,'cmbstyle','fullfast');
%freqnoise.crsspctrm = freqnoise.crsspctrm./(trace(freqnoise.crsspctrm)./size(freqnoise.crsspctrm,1));
%freqnoise2.crsspctrm = freqnoise2.crsspctrm./norm(freqnoise2.crsspctrm,'fro');

tmp   = ft_checkdata(freqnoise2, 'cmbstyle', 'fullfast');
denom = sqrt(norm(tmp.crsspctrm, 'fro'));
freqnoise2.fourierspctrm = freqnoise2.fourierspctrm./denom;

[a,b] = match_str(label, freqnoise.label);
lfpatch = lfpatch(a,:);
label   = label(a);

if numel(ndip)==1
  shuf = randperm(size(lfpatch,2),ndip);
else
  % treat ndip as a vector of indices for the active dipoles, the first two
  % are the interacting ones)
  shuf = ndip;
  ndip = numel(shuf);
end
lfpatch = lfpatch(:,shuf);
pos     = pos(shuf,:);
candidates = candidates(shuf);

% pos is assumed to be in mm, and leadfield.pos in m
pos = pos./1000;
ix  = zeros(size(pos,1),1);
iy  = ix;
insidevec = find(leadfield.inside);
for k = 1:size(pos,1)
  dpos = leadfield.pos-pos(k,:);
  d    = sqrt(sum(dpos.^2,2));
  [~,ix(k)] = min(d); 
  tmp = find(insidevec==ix(k));
  if ~isempty(tmp)
    iy(k) = tmp;
  else
    iy(k) = nan;
  end
end

if isempty(rho) && isempty(ampl)
  % new addition (May 2023) is to allow for non-scalar rho/amp/phs where
  % the vectors reflect elements of the lower triangular part of a matrix
  % of interacting dipoles: this allows for the modelling of more than 2
  % interacting dipoles. note that there is a constraint on how the phs
  % should behave, this is left at the discretion of the user, also the
  % vectors are not checked for their correct sizes at the moments
  ndip_int = (1+sqrt(1+8.*numel(amp)))./2;
  
  tri_low = tril(ones(ndip_int),-1)>0;
  tri_up  = triu(ones(ndip_int),1)>0;

  P = zeros(ndip_int);
  A = zeros(ndip_int);

  P(tri_low) = phs;
  P          = -P' + P;
  A(tri_low) = amp;
  A          = A' + A;

  rho      = eye(ndip);
  %rho(1,2) = amp.*exp(1i.*phs);
  %rho(2,1) = amp.*exp(-1i.*phs);
  %ampl     = [snr2 snr2 ones(1,ndip-2)-snr2]; %magnitude of source amplitude

  rho(1:ndip_int, 1:ndip_int) = eye(ndip_int) + A.*exp(1i.*P);
  ampl = [ones(1,ndip_int).*snr2 ones(1,ndip-ndip_int)-snr2];
% elseif ~isempty(rho) && ~isempty(ampl)
%   assert(all(size(rho)==ndip));
%   assert(all(size(ampl)==[1 ndip]));
else
  error('either rho and ampl should be both specified, or none of them');
end
freq           = freqnoise2;
%freq.crsspctrm = lfpatch*diag(ampl)*rho*diag(ampl)*lfpatch';
%freq.crsspctrm = freq.crsspctrm./(trace(freq.crsspctrm)./size(freq.crsspctrm,1));
%freq.crsspctrm = freq.crsspctrm./norm(freq.crsspctrm,'fro');

% due to the randomness, the output csd can deviate from the requested csd
for k = 1:500
  sx(k,1) = rng;
  [tmp, T] = mvnrnd(zeros(size(rho,1),1), diag(ampl)*rho*diag(ampl)', size(freqnoise2.fourierspctrm,1));
  %tmp = (randn(size(tmp))+1i.*randn(size(tmp)))*T;

  % rotate channel 1 and 2 a bit, to avoid that only channel 2 gets
  % non-zero imaginary part interactions with the rest, this is the
  % original implementation where only 2 connected dipoles are simulated
  tmp(:,1:ndip_int) = tmp(:,1:ndip_int).*exp(1i.*2.*pi.*rand(size(tmp,1),1));

  Cx = ctranspose(tmp)*tmp;
  Cx = Cx(1:ndip_int, 1:ndip_int);

  Aa(:,k) = diag(Cx);
  Cx = Cx./sqrt(abs(diag(Cx))*abs(diag(Cx))');
  rx(:,k) = abs(Cx(tri_low));
  ax(:,k) = angle(Cx(tri_low));
end

rx2 = abs(rx-amp(:))./std(rx-amp(:));
ax2 = abs(ax-phs(:))./std(ax-phs(:));
Ax2 = abs(mean(Aa(1:ndip_int,:))./size(freqnoise2.fourierspctrm,1)-snr2.^2);
Ax2 = Ax2./std(Ax2);

% get the point closest to the origin
d = sqrt(mean(rx2).^2+mean(ax2).^2+mean(Ax2).^2);
[mind, ixd] = min(d);
rng(sx(ixd));
tmp = mvnrnd(zeros(size(rho,1),1), diag(ampl)*rho*diag(ampl)', size(freqnoise2.fourierspctrm,1));
tmp(:,1:ndip_int) = tmp(:,1:ndip_int).*exp(1i.*2.*pi.*rand(size(tmp,1),1));
  
Ssim = (ctranspose(tmp)*tmp)./size(freqnoise2.fourierspctrm,1);
  
freq.fourierspctrm = (lfpatch*tmp')';
tmp   = ft_checkdata(freq, 'cmbstyle', 'fullfast');
denom = sqrt(norm(tmp.crsspctrm, 'fro'));
freq.fourierspctrm = freq.fourierspctrm./denom;
freq.fourierspctrm = snr.*freq.fourierspctrm + (1-snr).*freqnoise2.fourierspctrm;
freq      = ft_checkdata(freq,      'cmbstyle', 'fullfast');
freqnoise = ft_checkdata(freqnoise, 'cmbstyle', 'fullfast');

%freqnoise2 = ft_checkdata(freqnoise2, 'cmbstyle', 'fullfast');
%freq.crsspctrm = snr.*freq.crsspctrm + (1-snr).*freqnoise2.crsspctrm;

[a,b] = match_str(leadfield.label, label);
leadfield.leadfield(leadfield.inside) = cellrowselect(leadfield.leadfield(leadfield.inside),a);

params = struct('pos', pos, 'lf', lfpatch, 'indices', candidates, 'rho', rho, 'ampl', ampl, 'snr', snr, 'snr2', snr2, 'ix', ix, 'iy', iy, 'patchsigma', patchsigma, 'Ssim', Ssim);
