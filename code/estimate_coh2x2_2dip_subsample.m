function [output] = estimate_coh2x2_2dip_subsample(sourcemodel, freq, freqnoise, varargin)

% [COH] = ESTIMATE_COH2x2_2DIP(SOURCEMODEL, FREQ)
%
% estimate the coherence between all dipole-pairs.
% using a pairwise dipole spatial filter. This is an adjusted version of
% estimate_coh2x2_2dip, with a (slower) inner loop on sensor subsampling
% per refindx, which allows for the computation of more sophisticated
% measures of location and scatter (rather than mean and std, which seem to
% be non-robust).

lambda      = ft_getopt(varargin, 'lambda',      []);
refindx     = ft_getopt(varargin, 'refindx',     []);
outputflags = ft_getopt(varargin, 'outputflags', [1 1 1 1 1]);
doprewhiten = ~isempty(freqnoise);
subsample   = ft_getopt(varargin, 'subsample');
fixedori    = ft_getopt(varargin, 'fixedori', 'lambda1');
ori0        = ft_getopt(varargin, 'ori0', []); % reference orientation, against which the estimated ori is compared
getpow      = ft_getopt(varargin, 'outputpow', 0);
state       = ft_getopt(varargin, 'state', 'shuffle');
Tol         = ft_getopt(varargin, 'invTol', 1e-12);
fitmethod   = ft_getopt(varargin, 'fitmethod', 1); % 1: fit slope of a single line, 2: fit 2 lines
avgflag     = ft_getopt(varargin, 'avgflag', 1); % 1: compute estimate of mean and sem of coh based on trimming (yuen_dep), or 2: use normal mean and sem (faster) 
subsampletrials = ft_getopt(varargin, 'subsampletrials', 1);
lfproject2two_before = ft_getopt(varargin, 'lfproject2two_before', true);

switch fixedori
  case 'n.a.'
    oriflag = 0;
  case 'lambda1'
    oriflag = 1;
  case 'unitnoise'
    oriflag = 2;
  case 'pls2x2'
    oriflag = 3;
  case 'unitnoise4x4'
    oriflag = 4;
  case 'pls2x2b'
    oriflag = 5;
  otherwise
    error('unsupported option for ''fixedori''');
end

keepindividual = numel(refindx)==1;

if ~isempty(subsample)
  Nsens = subsample.N;
  nrand = subsample.nrand;
else
  Nsens = nan;
end
if ~isfinite(Nsens)
  nrand = 1;
end

% determine the output requested
assert(isequal(size(outputflags),[1 5]));
getcoh  = outputflags(1);
getcoh0 = outputflags(2); % fitted estimate of 'null'-coherence
getcohc = outputflags(3); % corrected coherence estimate, based on Wens' trick
getcoh1 = outputflags(4); % coherence obtained with a single dipole in a beamformer model
getcoh1_0 = outputflags(5); % fitted estimate of 'null'-coherence in a single dipole beamformer model

% compute cross-spectra
fprintf('computing cross-spectrum for data ...\n');
tmp = freq;
tmp = ft_checkdata(tmp, 'cmbstyle', 'fullfast');
C   = tmp.crsspctrm;
clear tmp;

% check whether trial subsampling can be performed: data needs fourier in input
if subsampletrials<1 && ~isfield(freq, 'fourierspctrm')
  error('subsampling of trials can only be performed with fourierspctrm in input');
end

if doprewhiten
  tmp = freqnoise;
  tmp = ft_checkdata(tmp, 'cmbstyle', 'fullfast');
  N   = tmp.crsspctrm;
  clear tmp;
else 
  N   = eye(size(C,1));
end

% setting some variables
inside  = sourcemodel.inside; if islogical(inside), inside = find(inside); end
ninside = numel(inside);
nchan   = numel(freq.label);

% reduce to two column leadfields, if needed
nori = size(sourcemodel.leadfield{inside(1)},2);
if nori==3 && lfproject2two_before
  fprintf('rotating leadfields into their 2-dimensional basis prior to subsampling\n');
  fprintf('and norm normalizing leadfields with the norm of the leadfield\n');
  lf = zeros(nchan, ninside*2);
  for k = 1:ninside
    tmplf   = sourcemodel.leadfield{inside(k)};
    [~,~,v] = svd(tmplf, 'econ');
    tmplf   = tmplf*v(:,1:2);
    lf(:,(k-1)*2+(1:2)) = tmplf./norm(tmplf);
  end
elseif nori==3
  lf = cat(2,sourcemodel.leadfield{inside});
elseif nori==2
  lf = cat(2,sourcemodel.leadfield{inside});
elseif nori==1
  assert(oriflag==0);
  lf = cat(2,sourcemodel.leadfield{inside});
end

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
dokappa = false;
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  lambda = ratio * trace(C)/size(C,1);
elseif ~isempty(lambda) && ischar(lambda) && lambda(end)=='#'
  dokappa = true;
  kappa  = sscanf(lambda,'%f#');
  lambda = 0;
elseif lambda>1
  % use lambda as 'kappa', i.e. do a truncated svd inverse
  dokappa = true;
  kappa  = lambda;
  lambda = 0;
elseif isempty(lambda)
  lambda = 0;
end

if isempty(refindx)
  refindx = 1:ninside;%100;
end
nref = numel(refindx);

% indices for sparse matrix that does the orientation projection
yproj = (1:2*ninside)';
xproj = reshape(repmat(1:ninside,[2 1]),       [], 1);
zproj = reshape([1:4:ninside*4;2:4:ninside*4], [], 1);

% preallocate memory for the 'workhorse variables'
lfClf   = zeros(2,2,ninside);
lfClfC  = zeros(2,2,ninside);
lfCClfC = zeros(2,2,ninside);

% preallocate memory for the output matrices
if getcoh,  coh  = zeros(ninside,nref,nrand) + 1i.*zeros(ninside,nref,nrand); end
if getcoh0
  coh0 = zeros(ninside,nref,nrand); 
  a     = zeros(nref,nrand);   
  r     = zeros(nref,nrand);
end
if getcohc, cohc = zeros(ninside,nref,nrand); end
if getcoh1, coh1 = zeros(ninside,nref,nrand); end
if getcoh1_0
  coh1_0 = zeros(ninside,nref,nrand);
  a1     = zeros(nref,nrand);
  r1     = zeros(nref,nrand);
end
if getpow
  sourcepow = zeros(ninside,nref,nrand); 
  powratio  = sourcepow;
end
if oriflag<=2 || oriflag==4
  ori = zeros(ninside,nrand);
end

rng(state);
nsub_ = zeros(nrand,1);
charcount = 0;
for m = 1:nrand
  cnt = 0;
  if subsampletrials<1
    % recompute the crosspectrum based on a smaller subset of trials
    ntrl = round(subsampletrials.*numel(freq.cumtapcnt));
    sel  = randperm(numel(freq.cumtapcnt), ntrl);

    tmpcfg = [];
    tmpcfg.trials = sel;
    tmp    = ft_selectdata(tmpcfg, freq);
    tmp    = ft_checkdata(tmp, 'cmbstyle', 'fullfast');
    C      = tmp.crsspctrm;
    
  end
  if isfinite(Nsens)
    % subsample C, N and leadfield
    if numel(Nsens)==1
      indx = sort(randperm(numel(freq.label),Nsens)); % keep it sorted!! -> subsampling of sensors
      fprintf(repmat('\b',1,charcount));
      charcount = fprintf('randomization %d/%d, number of channels is %d',m,nrand,numel(indx));
    elseif numel(Nsens)==2
      Ntmp = randi(Nsens(2)+1-Nsens(1),1)+Nsens(1)-1;
      indx = sort(randperm(numel(freq.label),Ntmp));
      fprintf(repmat('\b',1,charcount));
      charcount = fprintf('randomization %d/%d, number of channels is %d',m,nrand,numel(indx));
    elseif numel(Nsens)==nrand
      indx = sort(randperm(numel(freq.label),Nsens(k)));
    end
    nsub   = numel(indx);
    submat = sparse((1:nsub)',indx(:),ones(nsub,1),nsub,numel(freq.label));
    thisC  = submat*C*submat';
    thisN  = submat*N*submat';
    thislf = submat*lf;
    nsub_(m) = nsub;
  else
    thisN  = N;
    thisC  = C;
    thislf = lf;
  end

  if nori==3 && ~lfproject2two_before
    thislf3 = thislf;
    thislf  = zeros(size(thislf3,1), ninside*2);
    for k = 1:ninside
      tmplf   = thislf3(:,(k-1)*3+(1:3));
      [~,~,v] = svd(tmplf, 'econ');
      tmplf   = tmplf*v(:,1:2);
      thislf(:,(k-1)*2+(1:2)) = tmplf./norm(tmplf);
    end
  end

  if doprewhiten
    % compute the prewhitening matrix and apply this to C and leadfield
    [U,S,~] = svd(real(thisN));
    %Tol     = 1e-12;
    diagS   = diag(S);
    sel     = find(diagS>Tol.*diagS(1));
    P       = U(:,sel)*diag(1./sqrt(diag(S(sel,sel))))*U(:,sel)'; % prewhitening matrix
    
    thisC   = P*thisC*P';
    thislf  = P*thislf;
  else
    % don't prewhiten
  end
  
  rC    = real(thisC);
  rCreg = rC + eye(size(thisC,1))*lambda;
  
  if ~dokappa
    lfC     = rCreg\thislf; % lf'*inv(rCreg);
    if oriflag>=2, lfCsq = (rCreg^2)\thislf; end
  else
    fprintf('using kappa truncated svd for matrix inversion\n');
    [u,s,v] = svd(real(thisC));
    s = diag(s);
    if kappa<1
      thiskappa = round(numel(s)*kappa);
    else
      thiskappa = min(kappa, numel(s));
    end
    s = diag(1./s(1:thiskappa));
    invC = v(:,1:thiskappa)*s*u(:,1:thiskappa)';
    lfC  = invC*thislf;
    if oriflag>=2, lfCsq = (invC^2)*thislf; end
  end
  lfCC = thisC*lfC;
  
  if oriflag==0
    % nothing needs to be done for the orientation computation, already fixed
  else
    lfClf        = zeros(2,2,ninside);
    lfClf(1,1,:) = sum(lfC(:,1:2:end).*thislf(:,1:2:end));
    lfClf(2,1,:) = sum(lfC(:,1:2:end).*thislf(:,2:2:end));
    lfClf(1,2,:) = lfClf(2,1,:);
    lfClf(2,2,:) = sum(lfC(:,2:2:end).*thislf(:,2:2:end));

    lfCClfC(1,1,:) = abs(sum(lfCC(:,1:2:end).*lfC(:,1:2:end)));
    lfCClfC(2,1,:) = sum(lfCC(:,1:2:end).*lfC(:,2:2:end));
    lfCClfC(1,2,:) = conj(lfCClfC(2,1,:));
    lfCClfC(2,2,:) = abs(sum(lfCC(:,2:2:end).*lfC(:,2:2:end)));
  end

  if oriflag==1
    % compute fixedori based on max pow per single dipole
    pow   = sandwich2x2(inv2x2(lfClf), lfCClfC);
    [e,d] = eig2x2(real(pow));
    ori(:,m) = reshape(e(1,1,:),[],1)+1i.*reshape(e(2,1,:),[],1);
    if ~isempty(ori0)
      ori(:,m) = alignori(ori(:,m), ori0);
    end

  elseif oriflag==2
    % compute fixedori based on unitnoise gain
    lfCsqlf = zeros(size(lfClf));
    lfCsqlf(1,1,:) = sum(lfCsq(:,1:2:end).*thislf(:,1:2:end));
    lfCsqlf(2,1,:) = sum(lfCsq(:,1:2:end).*thislf(:,2:2:end));
    lfCsqlf(1,2,:) = lfCsqlf(2,1,:);
    lfCsqlf(2,2,:) = sum(lfCsq(:,2:2:end).*thislf(:,2:2:end));
    
    [e,d] = eig2x2(mtimes2x2(inv2x2(lfCsqlf),lfClf));
    ori(:,m) = reshape(e(1,1,:),[],1)+1i.*reshape(e(2,1,:),[],1);
  elseif oriflag>=3
    % this requires everything in '4x4'-space
  end
  
  if oriflag==0

    % project the lf-like matrices onto the fixedori
    [thislf, n] = normc(thislf);
    lfC         = (lfC)*spdiags(1./n(:),0,ninside,ninside);
    lfCC        = (lfCC)*spdiags(1./n(:),0,ninside,ninside);
  
    tmp1 = ( lfC(:,refindx(:))'*thislf)'; % if refindx is sparse, this is much faster and memory efficient than all-to-all
    tmp2 = (lfCC(:,refindx(:))'*lfC)';
    tmp3 = ( lfC(:,refindx(:))'*lfC)';
   
  elseif oriflag<3
    E     = sparse(yproj,xproj,e(zproj));

    % project the lf-like matrices onto the fixedori
    [thislf, n] = normc(thislf*E);
    lfC         = (lfC*E)*spdiags(1./n(:),0,ninside,ninside);
    lfCC        = (lfCC*E)*spdiags(1./n(:),0,ninside,ninside);
  
    tmp1 = ( lfC(:,refindx(:))'*thislf)'; % if refindx is sparse, this is much faster and memory efficient than all-to-all
    tmp2 = (lfCC(:,refindx(:))'*lfC)';
    tmp3 = ( lfC(:,refindx(:))'*lfC)';
  else
    % dipoles have 2 orientations
    refindx_ = [2.*refindx(:)-1 2.*refindx(:)]';
    tmp1 = ( lfC(:,refindx_(:))'*thislf)'; % if refindx is sparse, this is much faster and memory efficient than all-to-all
    tmp2 = (lfCC(:,refindx_(:))'*lfC)';
    tmp3 = ( lfC(:,refindx_(:))'*lfC)';
    tmp4 = (thislf(:,refindx_(:))'*thislf)';
    if oriflag==4
      tmp5 = (lfCsq(:,refindx_(:))'*thislf)';
    end
  end
  
  for k = refindx(:)'
    if oriflag<3
      cnt = cnt+1;
      % compute the lfClf etc. matrices, but now per fixedori dipole pair.
      if cnt==1
        % this is independent of the refindx
        lfClf(1,1,:)   = sum(lfC.*thislf);
        lfCClfC(1,1,:) = sum(lfCC.*conj(lfC));
        lfClfC(1,1,:)  = sum(lfC.*conj(lfC));
      end
      
      %lfClf(2,1,:) = sum(bsxfun(@times, lfC(:,k), thislf));
      lfClf(2,1,:) = tmp1(:,cnt).';
      lfClf(1,2,:) = conj(lfClf(2,1,:));
      lfClf(2,2,:) = sum(lfC(:,k).*thislf(:,k));
      denom        = inv2x2(lfClf);
      denom(:,:,k) = 1./lfClf(1,1,k);
      
      %lfCClfC(2,1,:) = sum(bsxfun(@times, lfCC(:,k), conj(lfC)));
      lfCClfC(2,1,:) = tmp2(:,cnt).';
      lfCClfC(1,2,:) = conj(lfCClfC(2,1,:));
      lfCClfC(2,2,:) = sum(lfCC(:,k).*conj(lfC(:,k)));
      numer          = lfCClfC;
      
      %lfClfC(2,1,:) = sum(bsxfun(@times, lfC(:,k), conj(lfC)));
      lfClfC(2,1,:) = tmp3(:,cnt).';
      lfClfC(1,2,:) = conj(lfClfC(2,1,:));
      lfClfC(2,2,:) = sum(lfC(:,k).*conj(lfC(:,k)));
      
      if getcohc || getcoh1 || getcoh1_0 || getpow
        % add (lf'*C*lf)^-1 to denom, as if it were a single dipole scan
        denom(3,1,:) = 1./lfClf(1,1,:);
        denom(4,2,:) = denom(3,1,k);
      end
      
      sel   = [1:k-1 k+1:ninside];
      csd   = sandwichMx2(denom, numer);
      
      pp          = sqrt(abs(csd(1,1,:)).*abs(csd(2,2,:)));
      coh(:,cnt,m)  = shiftdim(csd(1,2,:)./pp);
      coh(k,cnt,m)  = 1;
      
      if getpow
        sourcepow(:,cnt,m) = abs(csd(3,3,:));
      end
      
      if getcoh0
        csd0 = sandwich2x2(denom(1:2,1:2,:), lfClfC); %filter correlation, assuming diagonal noise
        
        % estimate 'scaling' parameter for null-csd
        if fitmethod==1
          [a(cnt,m),r(cnt,m)] = fit1line(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(1,2,sel)),[3 1 2]));
          %[a(sel,cnt,m),r(cnt,m)] = fitslope(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(1,2,sel)),[3 1 2]));
          coh0(sel,cnt, m) = (csd0(1,2,sel).*a(cnt,m))./pp(sel);
          %coh0(sel,cnt, m) = (shiftdim(csd0(1,2,sel)).*a(sel,cnt,m))./shiftdim(pp(sel));
        elseif fitmethod==2
          [a(cnt,m,1:2),r(cnt,m),idx] = fit2lines(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(1,2,sel)),[3 1 2]));
          coh0(sel( idx),cnt, m) = (csd0(1,2,sel( idx)).*a(cnt,m,1))./pp(sel( idx));
          coh0(sel(~idx),cnt, m) = (csd0(1,2,sel(~idx)).*a(cnt,m,2))./pp(sel(~idx));
        elseif fitmethod==3
          [a(cnt,m,1:2),r(cnt,m)] = fitcubic(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(1,2,sel)),[3 1 2]));
          coh0(sel, cnt, m) = (csd0(1,2,sel).*a(cnt,m,1) + (csd0(1,2,sel).^3).*a(cnt,m,2))./pp(sel);
        end
        coh0(k,  cnt, m) = 1;
      end
      
      if getcohc
        %csd = sandwich2x2(denom([1 4],:,:), numer);
        
        % this is equivalent to Wens' geometrical trick, where one of the
        % dipoles' leadfiedl is regressed out of the other one
        pp          = sqrt(abs(csd(1,1,:)).*abs(csd(4,4,:)));
        cohc(:,cnt,m) = csd(1,4,:)./pp;
        cohc(k,cnt,m) = 0;
      end
      
      if getcoh1
        %csd = sandwich2x2(denom([3 4],:,:), numer);
        
        % this is equivalent to a single dipole model in the beamformer
        pp = sqrt(abs(csd(3,3,:)).*abs(csd(4,4,:)));
        coh1(:,cnt,m) = csd(3,4,:)./pp;
      end
      
      if getcoh1_0
        csd0 = sandwich2x2(denom(3:4,1:2,:), lfClfC); %filter correlation, assuming diagonal noise
        
        % estimate 'scaling' parameter for null-csd
        if fitmethod==1
          [a1(cnt,m),r1(cnt,m)] = fit1line(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(3,4,sel)),[3 1 2]));
          %[a(sel,cnt,m),r(cnt,m)] = fitslope(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(1,2,sel)),[3 1 2]));
          coh1_0(sel,cnt, m) = (csd0(1,2,sel).*a1(cnt,m))./pp(sel);
          %coh0(sel,cnt, m) = (shiftdim(csd0(1,2,sel)).*a(sel,cnt,m))./shiftdim(pp(sel));
        elseif fitmethod==2
          [a1(cnt,m,1:2),r1(cnt,m),idx] = fit2lines(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(3,4,sel)),[3 1 2]));
          coh1_0(sel( idx),cnt, m) = (csd0(1,2,sel( idx)).*a1(cnt,m,1))./pp(sel( idx));
          coh1_0(sel(~idx),cnt, m) = (csd0(1,2,sel(~idx)).*a1(cnt,m,2))./pp(sel(~idx));
        elseif fitmethod==3
          [a1(cnt,m,1:2),r1(cnt,m)] = fitcubic(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(3,4,sel)),[3 1 2]));
          coh1_0(sel, cnt, m) = (csd0(1,2,sel).*a1(cnt,m,1) + (csd0(1,2,sel).^3).*a1(cnt,m,2))./pp(sel);
        end
        coh1_0(k,  cnt, m) = 1;
        
      end
    else
      cnt  = cnt+1;
      cnt_ = [2*cnt-1 2*cnt];
      k_   = [2*k-1   2*k];
      
      if cnt==1
        %lfClf2x2   = lfClf;
        %lfCClfC2x2 = lfCClfC;
        
        lfClf(4,4,end)   = 0;
        lfCClfC(4,4,end) = 0;
        
        lfClfC           = zeros(size(lfClf));
        lfClfC(1,1,:)    = sum(lfC(:,1:2:end).*conj(lfC(:,1:2:end)));
        lfClfC(2,1,:)    = sum(lfC(:,2:2:end).*conj(lfC(:,1:2:end)));
        lfClfC(1,2,:)    = sum(lfC(:,1:2:end).*conj(lfC(:,2:2:end)));
        lfClfC(2,2,:)    = sum(lfC(:,2:2:end).*conj(lfC(:,2:2:end)));
        
        lflf           = zeros(size(lfClf));
        lflf(1,1,:)    = sum(thislf(:,1:2:end).*conj(thislf(:,1:2:end)));
        lflf(2,1,:)    = sum(thislf(:,2:2:end).*conj(thislf(:,1:2:end)));
        lflf(1,2,:)    = sum(thislf(:,1:2:end).*conj(thislf(:,2:2:end)));
        lflf(2,2,:)    = sum(thislf(:,2:2:end).*conj(thislf(:,2:2:end)));
        
      end
      
      lfClf(3,1,:) = tmp1(1:2:end,cnt_(1)).';
      lfClf(4,1,:) = tmp1(1:2:end,cnt_(2)).';
      lfClf(3,2,:) = tmp1(2:2:end,cnt_(1)).';
      lfClf(4,2,:) = tmp1(2:2:end,cnt_(2)).';
      lfClf(1:2,3:4,:) = conj(permute(lfClf(3:4,1:2,:),[2 1 3]));
      lfClf(3,3,:) = sum(lfC(:,k_(1)).*thislf(:,k_(1)));
      lfClf(4,3,:) = sum(lfC(:,k_(2)).*thislf(:,k_(1)));
      lfClf(3,4,:) = sum(lfC(:,k_(1)).*thislf(:,k_(2)));
      lfClf(4,4,:) = sum(lfC(:,k_(2)).*thislf(:,k_(2)));
      denom        = inv4x4(lfClf);
      denom(:,:,k) = 1./lfClf(1,1,k);
      
      lfCClfC(3,1,:) = tmp2(1:2:end,cnt_(1)).';
      lfCClfC(4,1,:) = tmp2(1:2:end,cnt_(2)).';
      lfCClfC(3,2,:) = tmp2(2:2:end,cnt_(1)).';
      lfCClfC(4,2,:) = tmp2(2:2:end,cnt_(2)).';
      lfCClfC(1:2,3:4,:) = conj(permute(lfCClfC(3:4,1:2,:),[2 1 3]));
      lfCClfC(3,3,:) = abs(sum(lfCC(:,k_(1)).*conj(lfC(:,k_(1)))));
      lfCClfC(4,3,:) = sum(lfCC(:,k_(2)).*conj(lfC(:,k_(1))));
      lfCClfC(3,4,:) = conj(lfCClfC(4,3,:));
      %lfCClfC(3,4,:) = sum(lfCC(:,k_(1)).*conj(lfC(:,k_(2))));
      lfCClfC(4,4,:) = abs(sum(lfCC(:,k_(2)).*conj(lfC(:,k_(2)))));
      numer          = lfCClfC;
      
      lfClfC(3,1,:) = tmp3(1:2:end,cnt_(1)).';
      lfClfC(4,1,:) = tmp3(1:2:end,cnt_(2)).';
      lfClfC(3,2,:) = tmp3(2:2:end,cnt_(1)).';
      lfClfC(4,2,:) = tmp3(2:2:end,cnt_(2)).';
      lfClfC(1:2,3:4,:) = conj(permute(lfClfC(3:4,1:2,:),[2 1 3]));
      lfClfC(3,3,:) = sum(lfC(:,k_(1)).*conj(lfC(:,k_(1))));
      lfClfC(4,3,:) = sum(lfC(:,k_(1)).*conj(lfC(:,k_(2))));
      lfClfC(3,4,:) = sum(lfC(:,k_(2)).*conj(lfC(:,k_(1))));
      lfClfC(4,4,:) = sum(lfC(:,k_(2)).*conj(lfC(:,k_(2))));
      
      lflf(3,1,:) = tmp4(1:2:end,cnt_(1)).';
      lflf(4,1,:) = tmp4(1:2:end,cnt_(2)).';
      lflf(3,2,:) = tmp4(2:2:end,cnt_(1)).';
      lflf(4,2,:) = tmp4(2:2:end,cnt_(2)).';
      lflf(1:2,3:4,:) = conj(permute(lflf(3:4,1:2,:),[2 1 3]));
      lflf(3,3,:) = sum(thislf(:,k_(1)).*conj(thislf(:,k_(1))));
      lflf(4,3,:) = sum(thislf(:,k_(1)).*conj(thislf(:,k_(2))));
      lflf(3,4,:) = sum(thislf(:,k_(2)).*conj(thislf(:,k_(1))));
      lflf(4,4,:) = sum(thislf(:,k_(2)).*conj(thislf(:,k_(2))));
      
      sel   = [1:k-1 k+1:ninside];
      csd   = sandwich4x4(denom, numer);
      
      if oriflag==3 || oriflag==5
        x1 = mtimes2x2(csd(1:2,3:4,:),csd(3:4,1:2,:));
        x2 = mtimes2x2(csd(3:4,1:2,:),csd(1:2,3:4,:));
      
        [e1,d1]=eig2x2(real(x1));
        [e2,d2]=eig2x2(real(x2));
        
        assert(all(d1(1,1,:)>d1(2,2,:)));
        assert(all(d2(1,1,:)>d2(2,2,:)));
        
        projmat = zeros(2,4,ninside);
        projmat(1,1:2,:) = permute(e1(1:2,1,:),[2 1 3]);
        if oriflag==5
          % it seems that the orientations are constrained in the right
          % halfcircle anyway, so no orientation flip seems needed
          e2avg = mean(e2(1:2,1,:),3);
          e2(1:2,1,:) = repmat(e2avg./norm(e2avg),1,ninside);
        end
        projmat(2,3:4,:) = permute(e2(1:2,1,:),[2 1 3]);

      elseif oriflag==4
        if cnt==1
          lfCsqlf          = zeros(size(lfClf));
          lfCsqlf(1,1,:)   = sum(lfCsq(:,1:2:end).*conj(thislf(:,1:2:end)));
          lfCsqlf(2,1,:)   = sum(lfCsq(:,2:2:end).*conj(thislf(:,1:2:end)));
          lfCsqlf(1,2,:)   = sum(lfCsq(:,1:2:end).*conj(thislf(:,2:2:end)));
          lfCsqlf(2,2,:)   = sum(lfCsq(:,2:2:end).*conj(thislf(:,2:2:end)));
        end
        lfCsqlf(3,1,:) = tmp5(1:2:end,cnt_(1)).';
        lfCsqlf(4,1,:) = tmp5(1:2:end,cnt_(2)).';
        lfCsqlf(3,2,:) = tmp5(2:2:end,cnt_(1)).';
        lfCsqlf(4,2,:) = tmp5(2:2:end,cnt_(2)).';
        lfCsqlf(1:2,3:4,:) = conj(permute(lfCsqlf(3:4,1:2,:),[2 1 3]));
        lfCsqlf(3,3,:) = sum(lfCsq(:,k_(1)).*conj(thislf(:,k_(1))));
        lfCsqlf(4,3,:) = sum(lfCsq(:,k_(1)).*conj(thislf(:,k_(2))));
        lfCsqlf(3,4,:) = sum(lfCsq(:,k_(2)).*conj(thislf(:,k_(1))));
        lfCsqlf(4,4,:) = sum(lfCsq(:,k_(2)).*conj(thislf(:,k_(2))));
        
        tmp = mtimes4x4(inv4x4(lfCsqlf),lfClf);
        
        [e,d] = eig4x4(permute(tmp,[3 1 2]), 1);
        
        projmat = zeros(2,4,ninside);
        projmat(1,1:2,:) = permute(e(:,1:2), [2 1]);
        projmat(2,3:4,:) = permute(e(:,3:4), [2 1]);
        %error('this is probably nonsense');

        ori(:,m) = reshape(projmat(2,3,:),[],1)+1i.*reshape(projmat(2,4,:),[],1);
        if ~isempty(ori0)
          ori(:,m) = alignori(ori(:,m), ori0);
        end

      end
      
      
      
%  % this chunk of code does a 'spinning' approach to screen through
%  % pairs of orientations to compute the coh and coh0 brute force
%       ang  = linspace(-pi/2,pi/2,25);
%       nang = numel(ang);
%       
%       % the e matrix here was based on the first eigenvalue computation
%       for the single dipole orientation estimation
%       [offset,rho] = cart2pol(squeeze(e(1,1,:)),squeeze(e(2,1,:)));
%       
%       tmpcoh  = zeros(ninside,nang,nang,nref);
%       tmpcoh0 = zeros(ninside,nang,nang,nref);
%       tmplfc  = zeros(ninside,nang,nang,nref);
%       tmppow  = zeros(ninside,nang,nang,nref);
%       for kk = 1:numel(ang)
%         [ix,iy]=pol2cart(offset(:)+ang(kk),1);
%         for mm = 1:numel(ang)
%           [ix2,iy2]=pol2cart(ang(mm)+offset(k),1);
%         
%           projmat(1,1,:) = ix;
%           projmat(1,2,:) = iy;
%           projmat(2,3,:) = ix2;
%           projmat(2,4,:) = iy2;
%           
%           lfClf_   = sandwichMx4(projmat, lfClf);
%           lfCClfC_ = sandwichMx4(projmat, lfCClfC);
% %           lfClfC_  = sandwichMx4(projmat, lfClfC);
% %           lflf_    = sandwichMx4(projmat, lflf);
% %           
%           denom = inv2x2(lfClf_);
%           csd   = sandwich2x2(denom, lfCClfC_);
%           %csd0  = sandwich2x2(denom, lfClfC_);
%           pp    = sqrt(abs(csd(1,1,:)).*abs(csd(2,2,:)));
%       
%           % estimate 'scaling' parameter for null-csd
% %           [a(cnt,m),r(cnt,m)] = fitslope(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(1,2,sel)),[3 1 2]));
% %           
% %           A(kk,mm,cnt,m) = a(cnt,m);
% %           
% %           coh0(sel,cnt, m) = (csd0(1,2,sel).*a(cnt,m))./pp(sel);
% %           coh0(k,  cnt, m) = 1;
% %           
% %           tmpcoh(:,kk,mm,cnt)  = squeeze( csd(1,2,  :))./squeeze(pp);
% %           tmpcoh0(:,kk,mm,cnt) = squeeze(coh0(:,cnt,m));
% %           tmplfc(:,kk,mm,cnt)  = squeeze(lflf_(1,2,:))./squeeze(sqrt(lflf_(1,1,:).*lflf_(2,2,:)));
% %           tmppow(:,kk,mm,cnt,m)  = squeeze(abs(csd(1,1,:)));
% %           tmppow2(:,kk,mm,cnt,m) = squeeze(abs(csd(2,2,:)));
%           tmppow(:,kk,mm,cnt) = squeeze(pp);
%         end
%       end
%       
%       [maxval,ix]  = max(reshape(abs(tmppow),ninside,[]),[],2);
%       [i1,i2] = ind2sub([nang nang],ix);
%       offsetk = ang(i1);
%       offsetm = ang(i2);
%       
%       [ix,iy]   = pol2cart(offsetk(:)+offset(:),1);
%       [ix2,iy2] = pol2cart(offsetm(:)+offset(k),1);
%       
%       projmat(1,1,:) = ix;
%       projmat(1,2,:) = iy;
%       projmat(2,3,:) = ix2;
%       projmat(2,4,:) = iy2;
%       


      lfClf_   = sandwichMx4(projmat, lfClf);
      lfCClfC_ = sandwichMx4(projmat, lfCClfC);
      lfClfC_  = sandwichMx4(projmat, lfClfC);

      denom = inv2x2(lfClf_);
      csd   = sandwich2x2(denom, lfCClfC_);
      csd0  = sandwich2x2(denom, lfClfC_);
       
      pp            = sqrt(abs(csd(1,1,:)).*abs(csd(2,2,:)));
      coh(:,cnt,m)  = squeeze(csd(1,2,:)./pp);
      coh(k,cnt,m)  = 1;
      
      % estimate 'scaling' parameter for null-csd
      [a(cnt,m),r(cnt,m)] = fitslope(permute(real(csd0(1,2,sel)),[3 1 2]),permute(real(csd(1,2,sel)),[3 1 2]));
       
      coh0(sel,cnt, m) = squeeze((csd0(1,2,sel).*a(cnt,m))./pp(sel));
      coh0(k,  cnt, m) = 1;
      
      
      
        
    end % oriflag
    
    if getpow
      sourcepow(:,cnt,m) = shiftdim(abs(csd(1,1,:)));
      powratio(:,cnt,m)  = shiftdim(abs(csd(2,2,:))./abs(csd(1,1,:)));
    end
  end % for k = refindx(:)'
end % for m = 1:nrand

% compute more robust quantities across subsamples
if nrand>1
  if getcoh0
    if ~keepindividual && avgflag==1
      [~, coh, coh0, dvar] = yuent_dep(abs(coh), abs(coh0), 0.2, 3, 0);
    elseif ~keepindividual && avgflag==2
      dvar = nanstd(abs(coh)-abs(coh0),[],3)./sqrt(nrand);
      coh  = mean(abs(coh),3);
      coh0 = mean(abs(coh0),3);
    else
      coh  = squeeze(coh);
      coh0 = squeeze(coh0);
    end
  end

  if getpow
    if ~keepindividual
      sourcepow = mean(sourcepow,3);
      powratio  = mean(powratio,3);
    else
      sourcepow = squeeze(sourcepow);
      powratio  = squeeze(powratio);
    end
  end

  if getcohc
    if ~keepindividual
      cohc = mean(cohc,3); 
    else
      cohc = squeeze(cohc);
    end
  end
  
  if getcoh1 && getcoh1_0
    if ~keepindividual && avgflag==1
      [~, coh1, coh1_0, dvar1] = yuent_dep(abs(coh1), abs(coh1_0), 0.2, 3, 0);
    elseif ~keepindividual && avgflag==2
      dvar1  = nanstd(abs(coh1)-abs(coh1_0),[],3)./sqrt(nrand);
      coh1   = mean(abs(coh1),3);
      coh1_0 = mean(abs(coh1_0),3);
    else
      coh1 = squeeze(coh1);
      coh1_0 = squeeze(coh1_0);
    end
    
  elseif getcoh1
    if ~keepindividual && avgflag==1
      coh1 = median(coh1,3);
    elseif ~keepindividual && avgflag==2
      coh1 = mean(coh1,3);
    else
      coh1 = squeeze(coh1);
    end
  elseif getcoh1_0
    if ~keepindividual && avgflag==1
      coh1_0 = median(coh1_0,3);
    elseif ~keepindividual && avgflag==2
      coh1_0 = mean(coh1_0,3);
    else
      coh1_0 = squeeze(coh1_0);
    end
  end
end
  
fprintf('creating output structure\n');
if getcoh,  output.coh   = single(coh);  clear coh;  end
if getcoh0, output.coh0  = single(coh0); clear coh0; end
if getcohc, output.cohc  = single(cohc); clear cohc; end
if getcoh1, output.coh1  = single(coh1); clear coh1; end
if getcoh1_0, output.coh1_0 = single(coh1_0); clear coh1_0_; end
if exist('dvar1', 'var'), output.dvar1 = single(dvar1); clear dvar1; end
if exist('dvar', 'var'),  output.dvar  = single(dvar);  clear dvar;  end
if exist('ori', 'var'),   output.ori   = single(ori); end
if getpow
  output.pow      = single(sourcepow); clear sourcepow; 
  output.powratio = single(powratio);  clear powratio;
end

if getcoh0
  output.a = a;
  output.r = r;
end
if getcoh1_0 
  output.a1 = a1;
  output.r1 = r1;
end
if exist('ori','var')
  output.ori = ori;
end
output.refindx = refindx;
output.nsub    = nsub_(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use linear regression to fit the 'slope', i.e. noise scaling needed to
% get the coh0 and coh aligned
function [a,r,m,n] = fitslope(x,y)

xorig = x;

x = abs(x);
y = abs(y);

%[t1, t2] = percthreshold(x, 0.5, 1, 0);
%t1  = prctile(x, 95);

if 1
  t1  = median(x);%prctile(x, 50);
  sel = x>t1.*10;%.*0.5;
  
  x  = x(sel);
  mx = mean(x);
  
  y  = y(sel);
  my = mean(y);
  
  %w = spdiags(x(:), 0, n, n);
  x = x-mx;
  y = y-my;
  %a = (x'*w*y)./(x'*w*x);
  xsq = x.^2;
  a = (xsq'*y)./(xsq'*x); % this is the same as the weighted regression above, with weights the values of x
  r = sum((y-a*x).^2)./(y'*y);
  
  if nargout>2
    m = a*xorig;
    n = sum(sel);
  end
else
  [srt, ix] = sort(x);
  
  krn = ones(999,1)./999;
  xm_(ix,1) = convn(x(ix(:)),krn,'same');
  ym_(ix,1) = convn(y(ix(:)),krn,'same');
  
  x = x-xm_;
  y = y-ym_;
  
  xdenom(ix,1) = convn(x(ix(:)).^2,krn,'same');
  xy(ix,1)     = convn(x(ix(:)).*y(ix(:)),krn,'same');
  a(:,1)       = xy./xdenom; 
  r            = sum((y-a.*x).^2)./(y'*y);
  
end

function c = alignori(a, b)

% alignori aligns the complex values in a to have a minimum angle with the
% complex values in b
c = a;
a = a./abs(a);
b = b./abs(b);

ang = real(a).*real(b) + imag(a).*imag(b);
c(ang<0) = -c(ang<0);