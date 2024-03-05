

ndip = 20;
snr  = [0.6 0.5];
snr2 = 0.8;
amp  = [0.7 0.8];%[0 0.2 0.3 0.4 0.5 0.6];%[0.2 0.4 0.6];
phs  = 2.*pi.*[1/17 2/17 4/17 8/17];
patchsigma = 0;
N    = 200;

state = rng;
%state = output.state;
for i0 = 1:100
for i1 = 1:numel(ndip)
  for i2 = 1:numel(snr)
    for i3 = 1:numel(snr2)
      for i4 = 1:numel(amp)
        for i5 = 1:numel(phs)
          for i6 = 1:numel(patchsigma)
            for i7 = 1:numel(N) 
              %rng(state);
              
              batchid = sprintf('job_%03d_%d_%d_%d_%d_%d_%d_%d',i0,i1,i2,i3,i4,i5,i6,i7);
              
              qsubfeval2('run_script', 'script2', {'index_to_run' i0}, ...
                {'ndip' ndip(i1)}, ...
                {'snr'  snr(i2)},  ...
                {'snr2' snr2(i3)}, ...
                {'amp'  amp(i4)},  ...
                {'phs'  phs(i5)},  ...
                {'patchsigma' patchsigma(i6)}, ...
                {'N',   N(i7,:)},  ...
                'memreq', 8*1024^3, 'timreq', 359*60, 'batchid', batchid);
               pause(0.1);
            end
          end
        end
      end
    end
  end
end
end