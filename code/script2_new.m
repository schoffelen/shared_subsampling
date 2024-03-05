
if ~exist('index_to_run', 'var'), error('index_to_run needs to be specified'); end
if ~exist('state',      'var'), state      = rng(1, 'twister'); end % state of random number generator
if ~exist('ndip',       'var'), ndip       = 200; end % number of active dipoles
if ~exist('snr',        'var'), snr        = 0.5; end % brain-to-sensor relative SNR (snr vs. (1-snr))
if ~exist('snr2',       'var'), snr2       = 0.8; end % interacting source vs additional source snr (snr2 vs. (1-snr))
if ~exist('amp',        'var'), amp        = 0.8; end % amplitude of coupling (i.e. magnitude of coherence
if ~exist('phs',        'var'), phs        = 2.*pi./17; end % phase difference of coupling
if ~exist('N',          'var'), N          = 100;  end % number of sensors per subsample (can be 2 element)
if ~exist('nrand',      'var'), nrand      = 100;  end % number of subsamples
if ~exist('patchsigma', 'var'), patchsigma = 0.4;  end % spatial std of patch (determining size of active cortical patch
if ~exist('saveflag',   'var'), saveflag   = true; end % save output to a file
if ~exist('prewhiten',  'var'), prewhiten  = true; end % flag for prewhitening 
if ~exist('fitmethod',  'var'), fitmethod  = 1; end
if ~exist('noise', 'var'),      noise      = 'emptyroom'; end % added June 2023 to allow for toggling between emptyroom and task-data based noise
if ~exist('invTol', 'var'),     invTol     = 1e-12; end % added June 2023 to allow for configureable invTol
if ~exist('lambda', 'var'),     lambda     = 0; end % added June 2023 to allow for configureable lambda
if ~exist('lambda_sub', 'var'), lambda_sub = 0; end % added June 2023 to allow for configureable lambda_sub

for index_to_runx = index_to_run(:)'
  
  paramsfname = fullfile('/project/3011020.09/jansch_sandbox/tmp/data/params', sprintf('ndip%03d_%03d', ndip, index_to_runx));
  
  load(paramsfname);
  load('leadfield');
  
  [freq, freqnoise, params2] = simulate_data(struct('ndip',params.indx,'snr',snr,'snr2',snr2,'amp',amp,'phs',phs,'patchsigma',patchsigma,'noise',noise));
  if ~prewhiten
    freqnoise = [];
  end
  
  opts           = [];
  opts.N         = N;
  opts.nrand     = nrand;
  opts.ix        = params2.ix(:)';
  opts.prewhiten = prewhiten;
  opts.state     = state;
  opts.batchsize = 300; %100;
  opts.avgflag   = 2;
  opts.fitmethod = fitmethod;
  opts.ndip_int  = (1+sqrt(1+8.*numel(amp)))./2;
  opts.invTol    = invTol;
  opts.lambda    = lambda;
  opts.lambda_sub = lambda_sub;
  [coh, C, Cc, Corig, Corig1, c, c0, hl] = reconstruct_coh_new(freq, freqnoise, leadfield, opts);
  
  output        = struct('hl',hl,'params2',params2,'params',params,'opts',opts);
  output.C      = C;
  output.Cc     = Cc;
  output.Corig  = Corig;
  output.Corig1 = Corig1;
  output.state  = state;
  output.opts   = opts;
  
  if saveflag
    if numel(amp) == 1
      A = num2str(round(100*amp(1)),'%0.3d');
      P = num2str(round(100*phs(1)),'%0.3d');
    else
      A = 'amp_';
      P = 'phs_';
      for k = 1:numel(amp)
        A = cat(2, A, num2str(round(100*amp(k)), '%0.3d'), '_');
        P = cat(2, P, num2str(round(100*phs(k)), '%0.3d'), '_');
      end
      A = A(1:end-1);
      P = P(1:end-1);
    end

    outputfilename = strrep(paramsfname,'params','results');
    outputfilename = sprintf('%s_%0.3d_%0.3d_%s_%s_%0.3d_%0.3d',outputfilename,round(100*snr),round(100*snr2),A,P,round(100*patchsigma),opts.N(1));
    if ~prewhiten
      outputfilename = [outputfilename '_noprewhiten'];
    end
    save(outputfilename,'output');
  end
  
end