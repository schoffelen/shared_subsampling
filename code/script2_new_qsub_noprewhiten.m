

ndip = 20;
snr  = [0.6 0.5];
snr2 = 0.7;%[0.5 0.8];
amp  = [0 0.2 0.3 0.4 0.5 0.6 0.7 0.8];%[0.2 0.4 0.6];
phs  = 2.*pi.*[0 1/17 2/17 4/17 8/17];
patchsigma = 0;
N    = [50 150];%100;

state = rng(1,'twister');
%state = output.state;

list = reshape(1:100,[10 10])';

for i0 = 1:size(list,1)
for i1 = 1:numel(ndip)
  for i2 = 1:numel(snr)
    for i3 = 1:numel(snr2)
      for i4 = 1:numel(amp)
        for i5 = 1:numel(phs)
          for i6 = 1:numel(patchsigma)
            for i7 = 1:size(N,1) 
              %rng(state);
              
              batchid = sprintf('job_%03d_%d_%d_%d_%d_%d_%d_%d',i0,i1,i2,i3,i4,i5,i6,i7);
              
              qsubfeval('run_script', 'script2_new', {'index_to_run' list(i0,:)}, ...
                {'ndip' ndip(i1)}, ...
                {'snr'  snr(i2)},  ...
                {'snr2' snr2(i3)}, ...
                {'amp'  amp(i4)},  ...
                {'phs'  phs(i5)},  ...
                {'patchsigma' patchsigma(i6)}, ...
                {'N',   N(i7,:)},  ...
                {'nrand', 100}, ...
                {'prewhiten', false}, ...
                {'state', state}, ...
                {'fitmethod', 1}, ...
                'memreq', 16*1024^3, 'timreq', 479*60, 'batchid', batchid);
            
            %pause(0.1);
            end
          end
        end
      end
    end
  end
end
end