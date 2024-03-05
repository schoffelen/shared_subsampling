%load hits20210730_noprewhiten.mat

% get the intersection of all approaches that were successful
% allHits = Hits&HitsC&HitsOrig;
% sumHits = sum(allHits(:,:),2);
% [~,ix]  = max(sumHits);

% executive decision:
% indx = 15 (ix);
% snr = 0.6;
% phs = 8/17 pi
% amp = 0.6

extra_figs = 0;  % should extra figures beyond the ones used for he paper be plotted?


[~, whoami] = unix('whoami');

if contains(whoami, 'jansch')
    addpath /project/3011020.09/jansch_sandbox/tmp/
    addpath /project/3011020.09/jansch_sandbox/tmp/data
elseif contains(whoami, 'briwes')
    project_settings;  % set paths and add FieldTrip to path
    addpath(paths.code)  % github fork
else
    error('Do not know user.')
end


%% Load cmocean for perceptually uniform colormaps

% Download cmocean here:
% https://se.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps

addpath(paths.cmocean)


%% Run simulations or load output of simulations

if contains(whoami, 'jansch')

    %% this specifies the set of parameters used for the simulation of the example
    % that is used for most of the illustration figures.

    ndip          = 20;
    index_to_runx = 15;
    paramsfname   = fullfile('/project/3011020.09/jansch_sandbox/tmp/data/params', sprintf('ndip%03d_%03d', ndip, index_to_runx));
    prewhiten     = false;
    snr           = 0.6;
    snr2          = 0.8;
    amp           = 0.5;
    phs           = pi.*(8/17);
    patchsigma    = 0;
    N             = [50 150];%100;
    state         = rng(1,'twister');

    %% this is the 8mm resolution grid + forward model
    load(paramsfname);
    load('leadfield');
    leadfield_lowres = leadfield;


    %% this is the 4mm resolution grid + forward model, for the illustration
    resolution = 4;
    leadfield  = makesourcemodel_3dgrid(4);
    insidevec  = double(leadfield.inside);
    insidevec(insidevec>0) = 1:sum(insidevec);

    %% this simulates the sensor level cross-spectral density matrix, a noise
    % matrix (emptyroom) and a params2 structure, with details of the simulation
    % (indices of active dipoles etc)
    [freq, freqnoise, params2] = simulate_data(struct('ndip',params.indx,'snr',snr,'snr2',snr2,'amp',amp,'phs',phs,'patchsigma',patchsigma));

    %% this simulates sensor level data with the same activated sources, but now
    % with the interaction strength set to 0
    [freq0, ~, ~] = simulate_data(struct('ndip',params.indx,'snr',snr,'snr2',snr2,'amp',0,'phs',phs,'patchsigma',patchsigma));

    %% this simulates sensor level data with the same activated sources, but now
    % with the interaction between dipoles 1 and 2 at a phase diff of 0
    [freqphs0, freqnoisephs0, params2phs0] = simulate_data(struct('ndip',params.indx,'snr',snr,'snr2',snr2,'amp',amp,'phs',0,'patchsigma',patchsigma));

    % if ~prewhiten
    %   freqnoise = [];
    % end

    %% this computes an index array into the 4mm grid of the closest to the active
    % dipoles grid positions
    clear d ix
    for k = 1:size(params2.pos,1)
        [d(1,k),ix(1,k)] = min(sqrt(sum((leadfield.pos - params2.pos(k,:)).^2,2)));
    end

    %% this adds a few indices (into the 4mm grid) of points that are (close to)
    % the connecting line between the active dipoles
    % add some points along the line connecting dipoles 1 and 2
    L(:,1) = linspace(params2.pos(1,1),params2.pos(2,1),27);
    L(:,2) = linspace(params2.pos(1,2),params2.pos(2,2),27);
    L(:,3) = linspace(params2.pos(1,3),params2.pos(2,3),27);
    L = L(2:end-1,:);
    for k = 1:size(L,1)
        [d(end+1),ix(end+1)] = min(sqrt(sum((leadfield.pos - L(k,:)).^2,2)));
    end

    %% specify the options for the function that performs the subsampling,
    % specifically using the grid points close to the active dipoles as seeds,
    % as well as the ones on the connecting line between the interacting dipoles
    opts           = [];
    opts.N         = N;
    opts.nrand     = 100;
    opts.ix        = params2.ix;
    opts.prewhiten = prewhiten;
    opts.state     = state; % use the same random seed
    opts.batchsize = 177; %100;
    opts.refindx   = insidevec(ix); % should be expressed as indices, relative to the 'inside' only

    %% subsampling
    [coh, C, Cc, Corig, Corig1, c, c0, hl] = reconstruct_coh_new(freq, [], leadfield, opts);

    %% subsampling with interaction strength = 0
    [coh0, Cx, Ccx, Corigx, Corig1x, c00, c000, hlx] = reconstruct_coh_new(freq0, [], leadfield, opts);

    %% subsampling with interaction phase difference = 0
    [cohphs0, Cphs0, Ccphs0, Corigphs0, Corig1phs0, cphs0, c0phs0, hlphs0] = reconstruct_coh_new(freqphs0, [], leadfield, opts);

    % this is for evaluation of the subsampling for the 1-dipole case, one of
    % the missing links for a complete exposition -> probably not to be used
    opts.outputflags = [0 0 0 1 1];
    [coh2, C2, Cc2, Corig2, Corig2_1, c2, c20, hl2] = reconstruct_coh_new(freq, [], leadfield, opts);

    %[cohw, C, Cc, Corig, Corig1, cw, cw0, hl] = reconstruct_coh_new(freq, freqnoise, leadfield, opts);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOME FIGURES

    load('ndip020_015.mat');

    %%% JM 20230113
    % save the whole workspace here, for Britta
    mkdir('/project/3011020.09/jansch/tmp/forbritta');
    cd /project/3011020.09/jansch/tmp/forbritta
    save workspace_examples_20230113

elseif contains(whoami, 'briwes')

    load(paths.data)

    % change directory for saving figures
    cd(paths.figures)

end

%% FIG 1A

% spatial configuration of the activated dipoles.
s = leadfield_lowres;

visualize_cluster(data(2,2,5,4).C{8}(1), s, ...
    'highlight',{data(2,2,5,4).params2.pos(1:2,:); data(2,2,5,4).params2.pos(3:20,:)},...
    'highlcolor', [255, 127, 0]/255, ...
    'whitebg', 1, ...
    'plotval', 2)
set(gcf, 'color', 'w');
exportgraphics(gcf, 'dipoles.png');

%% Figure 3

comm = {};
nfp = cellfun(@numel, data(2,2,5,4).Corig);
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,4).Corig, s,'highlight',data(2,2,5,4).params2.pos(1:2,:), ...
    'comment', comm, 'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_2dip.png');

%% Figure 3A

% all-to-all pairwise coherence
% Single dipole beamformer

comm = {};
nfp = cellfun(@numel, data(2,2,5,4).Corig1);
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,4).Corig1, s, ...
    'highlight', data(2,2,5,4).params2.pos(1:2,:), ...
    'comment', comm, ...
    'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_orig1.png');

%% Figure 3B

% all-to-all pairwise coherence
% Subsampling two dipole beamformer

comm = {};
nfp = cellfun(@numel, data(2,2,5,4).C);
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,4).C, s,'highlight',data(2,2,5,4).params2.pos(1:2,:), ...
    'comment', comm, 'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_subsample.png');

%% Figure S2

comm = {};
nfp = cellfun(@numel, data(2,2,5,4).Cc); 
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,4).Cc, s,'highlight',data(2,2,5,4).params2.pos(1:2,:), 'comment', comm, ...
    'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_imagcoh.png');

%% 

comm = {};
nfp = cellfun(@numel, data(2,2,5,1).Corig);
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,1).Corig, s,'highlight',data(2,2,5,4).params2.pos(1:2,:), ...
    'comment', comm, 'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_2dip_phs0.png');

%%

comm = {};
nfp = cellfun(@numel, data(2,2,5,1).Corig1);
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,1).Corig1, s,'highlight',data(2,2,5,4).params2.pos(1:2,:), ...
    'comment', comm, 'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_orig1_phs0.png');

%% Figure 4B

comm = {};
nfp = cellfun(@numel, data(2,2,5,1).C);
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,1).C, s,'highlight',data(2,2,5,4).params2.pos(1:2,:), ...
    'comment', comm, 'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_subsample_phs0.png');

%% Figure 4A

comm = {};
nfp = cellfun(@numel, data(2,2,5,1).Cc);
thr = 1-[95 99 99.5 99.9 99.95 99.99 99.995 99.999 99.9995]/100;
for k = 1:9
    comm{k}{1,1} = sprintf('threshold: %1.5g', thr(k));
    comm{k}{2,1} = sprintf('# of conn: %d',nfp(k));
    comm{k}      = char(comm{k});
end

visualize_cluster(data(2,2,5,1).Cc, s,'highlight',data(2,2,5,4).params2.pos(1:2,:), ...
    'comment', comm, 'whitebg', 1);
set(gcf, 'color', 'w');
set(gcf, 'position', [50 50 3*560 3*420]);
set(gca, 'position', [0.01 0.01 0.98 0.98]);
exportgraphics(gcf, 'blobs_imagcoh_phs0.png');

%% Setup for Figure 1

% show the coherence, the null coherence, and the coherence for freq0

[volidx(:,2), volidx(:,1), volidx(:,3)] = ind2sub(leadfield.dim,ix); % swap is intentional
in = find(leadfield.inside);
out = find(~leadfield.inside);
dum = zeros(leadfield.dim);


% the volplots in the subplots is a 4x5 arrangement for slices 1:2:end-1
slice_idx = 1:2:39;
dim = size(dum); 
for k = 1:size(volidx,1)
    selz(k,1) = nearest(slice_idx, volidx(k,3));
    subx(k,1) = (mod(selz(k,1)-1, 5))*dim(1)+volidx(k,2); % column id
    suby(k,1) = (ceil(selz(k,1)/5)-1)*dim(2)+volidx(k,1); % row id
end

% colormaps
cmap_pos = cmocean('thermal');
cmap_pos(1:10, :) = [];
cmap_diff = cmocean('balance');

% exporting the graphics:
kwargs = {'resolution', 600, ...
    'backgroundcolor', 'current'};

%% Figure 1B, 1E, 1F

% 1B - Estimate for true connectivity
figure;
dip_indx=1;
dum(in)=abs(coh.coh1(:,dip_indx));
dum(out) = NaN;  % make it possible to plot background as white
% subplot('position',[0.005 0.505 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ws', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ys', 'markersize', 7, 'linewidth', 1.5)
caxis([0 1]);
ft_colormap(cmap_pos);
set(gcf, 'color', [1 1 1])
colorbar
exportgraphics(gcf, 'coh_seed1_nullcoh_fig0.png', kwargs{:});


% Figure 1E - Estimate of null coherence
figure;
dum(in)=abs(coh.coh1_0(:,dip_indx));
% subplot('position',[0.505 0.505 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ws', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ys', 'markersize', 7, 'linewidth', 1.5)
caxis([0 1]);
ft_colormap(cmap_pos);
set(gcf, 'color', [1 1 1])
colorbar
exportgraphics(gcf, 'coh_seed1_nullcoh_fig1.png', kwargs{:});


% Figure 1F - Difference map B - E
figure;
dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh.coh1_0(:,dip_indx));
% subplot('position',[0.2525 0.005 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ks', 'markersize', 7, 'linewidth', 1.5)
caxis([-.3 .3]);
ft_colormap(cmap_diff);
set(gcf, 'color', [1 1 1])

colorbar
exportgraphics(gcf, 'coh_seed1_nullcoh_fig2.png', kwargs{:});


%% Figures for 0 phase null coherence

if(extra_figs)
    figure;  %#ok
    dip_indx=1;
    dum(in)=abs(cohphs0.coh1(:,dip_indx));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([0 1]);
    ft_colormap(cmap_pos);
    set(gcf, 'color', 'w');
    colorbar


    %     figure;
    dum(in)=abs(cohphs0.coh1_0(:,dip_indx));
    subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([0 1]);
    ft_colormap(cmap_pos);
    set(gcf, 'color', 'w');
    colorbar

    %     figure;
    dum(in)=abs(cohphs0.coh1(:,dip_indx))-abs(cohphs0.coh1_0(:,dip_indx));
    subplot('position',[0.2525 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    set(gcf, 'color', 'w');
    colorbar

    exportgraphics(gcf, 'coh_seed1_nullcoh_phs0.png');
end

%% Figures 1C and 1D

figure;
dip_indx=1;
dum(in)=abs(coh.coh1(:,dip_indx));
% subplot('position',[0.005 0.505 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ws', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ys', 'markersize', 7, 'linewidth', 1.5);
caxis([0 1]);
ft_colormap(cmap_pos);
set(gcf, 'color', [1, 1, 1])
colorbar
exportgraphics(gcf, 'coh_seed1_coh0_fig0.png', kwargs{:});

% Figure 1C
figure
dum(in)=abs(coh0.coh1(:,dip_indx));
% subplot('position',[0.505 0.505 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ws', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ys', 'markersize', 7, 'linewidth', 1.5);
caxis([0 1]);
ft_colormap(cmap_pos);
colorbar
set(gcf, 'color', [1, 1, 1])
exportgraphics(gcf, 'coh_seed1_coh0_0_fig1.png', kwargs{:});

% Figure 1D
figure
dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh0.coh1(:,dip_indx));
% subplot('position',[0.2525 0.005 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ks', 'markersize', 7, 'linewidth', 1.5)
caxis([-.3 .3]);
ft_colormap(cmap_diff);
set(gcf, 'color', 'w');
colorbar
set(gcf, 'color', [1, 1, 1])
exportgraphics(gcf, 'coh_seed1_coh0_fig2.png', kwargs{:});

%% Figures for 2 dipole BF null coherence

if(extra_figs)
    figure;  %#ok
    dip_indx=1;
    dum(in)=abs(coh.coh(:,dip_indx));
    % subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([0 1]);
    ft_colormap(cmap_pos);
    colorbar

    figure
    dum(in)=abs(coh.coh0(:,dip_indx));
    % subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;plot(subx(1),suby(1),'ks');plot(subx(2),suby(2),'ys')
    caxis([0 1]);
    ft_colormap(cmap_pos);
    colorbar

    figure
    dum(in)=abs(coh.coh(:,dip_indx))-abs(coh.coh0(:,dip_indx));
    % subplot('position',[0.2525 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    set(gcf, 'color', 'w');
    colorbar

    exportgraphics(gcf, 'coh_seed1_2dip_nullcoh.png');

end

%% Figures for 0 phase 2 dipole BF null coherence

if(extra_figs)

    figure;  %#ok
    dip_indx=1;
    dum(in)=abs(cohphs0.coh(:,dip_indx));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([0 1]);
    ft_colormap(cmap_diff);
    colorbar

    %     figure
    dum(in)=abs(cohphs0.coh0(:,dip_indx));
    subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([0 1]);
    ft_colormap(cmap_diff);
    colorbar

    %     figure
    dum(in)=abs(cohphs0.coh(:,dip_indx))-abs(cohphs0.coh0(:,dip_indx));
    subplot('position',[0.2525 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    set(gcf, 'color', 'w');
    colorbar

    exportgraphics(gcf, 'coh_seed1_2dip_nullcoh_phs0.png');
end

%% Figures for subsampling 2 dipole BF null coherence

if(extra_figs)
    figure;  %#ok
    dip_indx=1;
    dum(in)=abs(c(:,dip_indx));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([0 15]);
    ft_colormap(cmap_pos);


    dum(in)=abs(c0(:,dip_indx));
    subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([0 15]);
    ft_colormap(cmap_pos);


    dum(in)=abs(c(:,dip_indx))-abs(c0(:,dip_indx));
    subplot('position',[0.2525 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([-6 6]);
    ft_colormap(cmap_diff);
    set(gcf, 'color', 'w');
    exportgraphics(gcf, 'coh_seed1_2dip_nullcoh_subsample.png');

end

%% Figures for 0 phase subsampling 2 dipole BF null coherence

if(extra_figs)

    figure;  %#ok
    dip_indx=1;
    dum(in)=abs(cphs0(:,dip_indx));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([0 15]);
    ft_colormap(cmap_pos);

    dum(in)=abs(c0phs0(:,dip_indx));
    subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([0 15]);
    ft_colormap(cmap_pos);

    dum(in)=abs(cphs0(:,dip_indx))-abs(c0phs0(:,dip_indx));
    subplot('position',[0.2525 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ys')
    caxis([-6 6]);
    ft_colormap(cmap_diff);
    set(gcf, 'color', 'w');
    exportgraphics(gcf, 'coh_seed1_2dip_nullcoh_subsample_phs0.png');

end

%% Figures for imaginary coherence

if(extra_figs)

    figure;  %#ok
    dip_indx=1;
    dum(in)=abs(imag(coh.cohc(:,dip_indx)));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_pos);
    set(gcf, 'color', 'w');
    exportgraphics(gcf, 'coh_seed1_imagcoh.png');

    %     figure;
    dip_indx=1;
    dum(in)=abs(imag(cohphs0.cohc(:,dip_indx)));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    set(gcf, 'color', 'w');
    exportgraphics(gcf, 'coh_seed1_imagcoh_phs0.png');

end

%% Figures Supplementary: different seeds

% show some difference maps for other seeds

% FIG S1A
figure
dip_indx=2;
dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh.coh1_0(:,dip_indx));
% subplot('position',[0.005 0.505 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ks', 'markersize', 7, 'linewidth', 1.5);
caxis([-.3 .3]);
ft_colormap(cmap_diff);
colorbar
set(gcf, 'color', [1, 1, 1])
exportgraphics(gcf, 'coh_seed24828_nullcoh_fig0.png', kwargs{:});


% FIG S1B
figure
dip_indx=4;
dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh.coh1_0(:,dip_indx));
% subplot('position',[0.505 0.505 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(4),suby(4),'ko', 'markerfacecolor', 'y', 'markersize', 6);
caxis([-.3 .3]);
ft_colormap(cmap_diff);
colorbar
set(gcf, 'color', [1, 1, 1])
exportgraphics(gcf, 'coh_seed24828_nullcoh_fig1.png', kwargs{:});


% FIG S1C
figure
dip_indx=8;
dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh.coh1_0(:,dip_indx));
% subplot('position',[0.005 0.005 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ks', 'linewidth', 1.5);
plot(subx(2),suby(2),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(8),suby(8),'ko', 'markerfacecolor', 'y', 'markersize', 6);
caxis([-.3 .3]);
ft_colormap(cmap_diff);
colorbar
set(gcf, 'color', [1, 1, 1])
exportgraphics(gcf, 'coh_seed24828_nullcoh_fig2.png', kwargs{:});


% FIG S1D
figure
dip_indx=28;
dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh.coh1_0(:,dip_indx));
% subplot('position',[0.505 0.005 0.49 0.49]);
volplot(dum(:,:,1:2:end-1),'montage');
hold on;
plot(subx(1),suby(1),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(2),suby(2),'ks', 'markersize', 7, 'linewidth', 1.5);
plot(subx(28),suby(28),'ko', 'markerfacecolor', 'y', 'markersize', 6);
caxis([-.3 .3]);
ft_colormap(cmap_diff);
colorbar
set(gcf, 'color', [1, 1, 1])
exportgraphics(gcf, 'coh_seed24828_nullcoh_fig3.png', kwargs{:});

%% Figures for extra seeds, 2 dipole BF

if(extra_figs)

    figure;  %#ok
    dip_indx=2;
    dum(in)=abs(coh.coh(:,dip_indx))-abs(coh.coh0(:,dip_indx));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    exportgraphics(gcf, 'coh_seed24828_2dip_nullcoh_fig0.png');


    dip_indx=4;
    dum(in)=abs(coh.coh(:,dip_indx))-abs(coh.coh0(:,dip_indx));
    subplot('position',[0.505 0.505 0.49 0.49]);v
    olplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(4),suby(4),'kp');
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    exportgraphics(gcf, 'coh_seed24828_2dip_nullcoh_fig1.png');


    dip_indx=8;
    dum(in)=abs(coh.coh(:,dip_indx))-abs(coh.coh0(:,dip_indx));
    subplot('position',[0.005 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(8),suby(8),'kp');
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    exportgraphics(gcf, 'coh_seed24828_2dip_nullcoh_fig2.png');


    dip_indx=28;
    dum(in)=abs(coh.coh(:,dip_indx))-abs(coh.coh0(:,dip_indx));
    subplot('position',[0.505 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(28),suby(28),'kp');
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    exportgraphics(gcf, 'coh_seed24828_2dip_nullcoh_fig3.png');

end

%% Figures for extra seeds

if(extra_figs)

    figure;  %#ok
    dip_indx=2;
    dum(in)=abs(imag(coh.cohc(:,dip_indx)));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_pos);

    dip_indx=4;
    dum(in)=abs(imag(coh.cohc(:,dip_indx)));
    subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(4),suby(4),'kp');
    caxis([-.3 .3]);
    ft_colormap(cmap_pos);

    dip_indx=8;
    dum(in)=abs(imag(coh.cohc(:,dip_indx)));
    subplot('position',[0.005 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(8),suby(8),'kp');
    caxis([-.3 .3]); ft_colormap(cmap_pos);

    dip_indx=28;
    dum(in)=abs(imag(coh.cohc(:,dip_indx)));
    subplot('position',[0.505 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(28),suby(28),'kp');
    caxis([-.3 .3]); ft_colormap(cmap_pos);
    exportgraphics(gcf, 'coh_seed24828_2dip_imagcoh.png');

end

%% Figures for  extra seed, 2 dipole - subsampling BF

if(extra_figs)

    figure;  %#ok
    dip_indx=2;
    dum(in)=abs(c(:,dip_indx))-abs(c0(:,dip_indx));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-6 6]);
    ft_colormap(cmap_pos);

    dip_indx=4;
    dum(in)=abs(c(:,dip_indx))-abs(c0(:,dip_indx));
    subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(4),suby(4),'kp');
    caxis([-6 6]);
    ft_colormap(cmap_pos);

    dip_indx=8;
    dum(in)=abs(c(:,dip_indx))-abs(c0(:,dip_indx));
    subplot('position',[0.005 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(8),suby(8),'kp');
    caxis([-6 6]);
    ft_colormap(cmap_pos);

    dip_indx=28;
    dum(in)=abs(c(:,dip_indx))-abs(c0(:,dip_indx));
    subplot('position',[0.505 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(28),suby(28),'kp');
    caxis([-6 6]);
    ft_colormap(cmap_pos);
    exportgraphics(gcf, 'coh_seed24828_2dip_nullcoh_subsample.png');
end

%%  Figures for extra seeds, no coherence

if(extra_figs)

    figure;  %#ok
    dip_indx=2;
    dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh0.coh1(:,dip_indx));
    subplot('position',[0.005 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks')
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);

    dip_indx=4;
    dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh0.coh1(:,dip_indx));
    subplot('position',[0.505 0.505 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(4),suby(4),'kp');
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);

    dip_indx=8;
    dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh0.coh1(:,dip_indx));
    subplot('position',[0.005 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(8),suby(8),'kp');
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);

    dip_indx=28;
    dum(in)=abs(coh.coh1(:,dip_indx))-abs(coh0.coh1(:,dip_indx));
    subplot('position',[0.505 0.005 0.49 0.49]);
    volplot(dum(:,:,1:2:end-1),'montage');
    hold on;
    plot(subx(1),suby(1),'ks');
    plot(subx(2),suby(2),'ks');
    plot(subx(28),suby(28),'kp');
    caxis([-.3 .3]);
    ft_colormap(cmap_diff);
    exportgraphics(gcf, 'coh_seed24828_coh0.png');

end

%% Figure 2 - Scatter plots for approaches (Nike plots)

% Those colors are colorblind safe (Prot./Deut./Trit.)
colors.single = [  0, 171, 142] ./ 255;  % #00AB8E
colors.twodip = [ 1, 114, 189] ./ 255;  % #0172BD
colors.subsam = [148,  23,  81] ./ 255; % #941751
colors.imag = [191, 130, 20] ./ 255;  % #BF8214

% scatter plots for a few dipoles
if(extra_figs)
    D = [1 2 4 8 28 32];  %#ok
    m_cases = 1:8;
else
    D = [1 2 32];
    m_cases = [1, 3, 4];
end

mrk = [insidevec(ix(2)) insidevec(ix(1)) 0 0 0 0];
for k = 1:numel(D)
    dip_indx = D(k);
    for m = m_cases
        figure;hold on;
        switch m
            case 1
                dat0 = abs(coh.coh1_0(1:2:end,dip_indx));
                dat1 = abs(coh.coh1(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(coh.coh1_0(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(coh.coh1(mrk(k), dip_indx)); end
                xlab = '''null'' coherence';
                col = colors.single;
                txt  = sprintf('single dipole beamformer, refdip %d', dip_indx);
                fname = sprintf('scatter_seed%d_singledipole_nullcoh.png', dip_indx);
            case 2
                dat0 = abs(coh0.coh1(1:2:end,dip_indx));
                dat1 = abs(coh.coh1(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(coh0.coh1(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(coh.coh1(mrk(k), dip_indx)); end
                xlab = 'coherence in null-data';
                txt  = sprintf('single dipole beamformer, refdip %d', dip_indx);
                fname = sprintf('scatter_seed%d_singledipole_coh0.png', dip_indx);
            case 3
                dat0 = abs(coh.coh0(1:2:end,dip_indx));
                dat1 = abs(coh.coh(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(coh.coh0(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(coh.coh(mrk(k), dip_indx));  end
                xlab = '''null'' coherence';
                col = colors.twodip;
                txt  = sprintf('two dipole beamformer, refdip %d', dip_indx);
                fname = sprintf('scatter_seed%d_doubledipole_nullcoh.png', dip_indx);
            case 4
                dat0 = abs(c0(1:2:end,dip_indx));
                dat1 = abs(c(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(c0(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(c(mrk(k), dip_indx)); end
                xlab = '''null'' coherence, Z-scored';
                col = colors.subsam;
                txt  = sprintf('subsampling two dipole beamformer, refdip %d', dip_indx);
                fname = sprintf('scatter_seed%d_doubledipole_nullcoh_subsample.png', dip_indx);
            case 5
                dat0 = abs(c20(1:2:end,dip_indx));
                dat1 = abs(c2(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(c20(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(c2(mrk(k), dip_indx)); end
                xlab = '''null'' coherence, Z-scored';
                txt  = sprintf('subsampling single dipole beamformer, refdip %d', dip_indx);
                fname = sprintf('scatter_seed%d_singledipole_nullcoh_subsample.png', dip_indx);
            case 6
                dat0 = abs(cohphs0.coh1_0(1:2:end,dip_indx));
                dat1 = abs(cohphs0.coh1(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(cohphs0.coh1_0(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(cohphs0.coh1(mrk(k), dip_indx)); end
                xlab = '''null'' coherence';
                txt  = sprintf('single dipole beamformer, refdip %d, phs0', dip_indx);
                fname = sprintf('scatter_seed%d_singledipole_nullcoh_phs0.png', dip_indx);
            case 7
                dat0 = abs(cohphs0.coh0(1:2:end,dip_indx));
                dat1 = abs(cohphs0.coh(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(cohphs0.coh0(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(cohphs0.coh(mrk(k), dip_indx)); end
                xlab = '''null'' coherence';
                txt  = sprintf('two dipole beamformer, refdip %d, phs0', dip_indx);
                fname = sprintf('scatter_seed%d_doubledipole_nullcoh_phs0.png', dip_indx);
            case 8
                dat0 = abs(c0phs0(1:2:end,dip_indx));
                dat1 = abs(cphs0(1:2:end,dip_indx));
                if mrk(k), ref0 = abs(c0phs0(mrk(k), dip_indx)); end
                if mrk(k), ref1 = abs(cphs0(mrk(k), dip_indx)); end
                xlab = '''null'' coherence';
                txt  = sprintf('subsampling two dipole beamformer, refdip %d, phs0', dip_indx);
                fname = sprintf('scatter_seed%d_doubledipole_nullcoh_phs0.png', dip_indx);
        end

        if isempty(col)
            col = [1, 114, 189 ]./ 255;
        end

        plot(dat0, dat1,'o', 'color', col)
        if mrk(k)
            plot(ref0, ref1, 'ks', 'markerfacecolor', [255 255 51] ./ 255, ...
                'markersize', 10);
        end
        axis equal;
        axis square
        set(gcf,'color','white')
        if m==4
            ylabel('coherence, z-scored');
        else
            ylabel('coherence');
        end
        if m==4 || m==5 || m==8
            plot([0 18], [0 18], 'k', 'linewidth', 1.5);
            set(gca,'tickdir','out','xtick',[0 4 8 12 16],'ytick',[0 4 8 12 16]);
            axis([0 18 0 18]);
        else
            plot([0 1], [0 1],'k', 'linewidth', 1.5);
            set(gca,'tickdir','out','xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0 0.2 0.4 0.6 0.8 1]);
            axis([0 1 0 1]);
        end
        xlabel(xlab)
        title(txt);
        if m==1
            abc = axis;
        elseif m==2
            axis(abc)
        end

        exportgraphics(gcf, fname, kwargs{:});
    end
end

%% Setup for figures 5 and 6

if contains(whoami, 'jansch')
    datadir = '/project/3011020.09/jansch/tmp/results_noprewhiten_20220615_020';
    %datadir = '/project/3011020.09/jansch/tmp/results_prewhiten_20220627_020';
elseif contains(whoami, 'briwes')
    datadir = paths.bargraph_data;
end

cd(datadir);
d = dir('*summary.mat');
for k = 1:numel(d)
    load(d(k).name);
    allC(k) = C; % subsampling
    allCc(k) = Cc; % imag coh
    allCorig(k) = Corig; % double dipole no subsampling
    allCorig1(k) = Corig1; % single dipole no subsampling
end
phs     = round(180*uvals{4}./(pi*100));
sel     = 1:100;

if contains(whoami, 'briwes')
    % change directory back for saving figures
    cd(paths.figures)
end

snrs = 1:2;  

% colors - those are colorblind safe (Prot./Deut./Trit.)
bar_colors = [
    0,   0, 139  % darkblue  #00008B
    80, 119, 189  % blue  #5077BD
    101, 204, 169  % aquarmarine  #65CCA9
    230, 160,  29  % gold  #E6A91D
    165, 104,  40, % brown  #A56828
    ] ./ 255;

%% Figure 5 - bar graphs

for i1 = snrs  % SNR
    for i2 = 1:3  % amplitude
        figure('position',[10 10 1120 840]);
        dx = permute(cat(6, allCorig1.D), [6 1 2 3 4 5]); % distance metric, goes to 0
        % when the simulated dipoles are within the 'blobs'

        %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;

        % the sum is across the 100 different dipole configurations
        % the max is across the set of thresholds tested, so the number
        % returned reflects the hit rate for the best threshold.
        Hx = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
        subplot(2,2,1);
        b = bar(uvals{3}./100, Hx);
        update_bar_colors(b, bar_colors)
        ylim([0 1]); ylabel('hit rate'); xlabel('coherence strength');
        legend(cellstr(num2str(phs)),'location','NorthWest');
        title('Coherence, single dipole BF')

        dx = permute(cat(6, allCorig.D), [6 1 2 3 4 5]);
        %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
        Hx = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
        subplot(2,2,2);
        b = bar(uvals{3}./100, Hx);
        update_bar_colors(b, bar_colors)
        ylim([0 1]); ylabel('hit rate'); xlabel('coherence strength');
        title('Coherence, two-dipole BF')

        dx = permute(cat(6, allCc.D), [6 1 2 3 4 5]);
        %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
        Hx = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
        subplot(2,2,3);
        b = bar(uvals{3}./100, Hx);
        update_bar_colors(b, bar_colors)
        ylim([0 1]); ylabel('hit rate'); xlabel('coherence strength');
        title('imag(coh), single dipole BF')

        dx = permute(cat(6, allC.D), [6 1 2 3 4 5]);
        %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
        Hx = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
        subplot(2,2,4);
        b = bar(uvals{3}./100, Hx);
        update_bar_colors(b, bar_colors)
        ylim([0 1]); ylabel('hit rate'); xlabel('coherence strength');
        title('Coherence, two-dipole subsampling BF')

        set(gcf, 'color', 'w')
        exportgraphics(gcf, sprintf('detection_rate_%d_%d.pdf', uvals{1}(i1), uvals{2}(i2)));
    end
end


if (extra_figs)
    for i1 = snrs  %#ok
        for i2 = 1:3
            figure('position',[10 10 1120 840]);
            dx = permute(cat(6, allCorig1.D), [6 1 2 3 4 5]); % distance metric, goes to 0
            % when the simulated dipoles are within the 'blobs'

            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;

            % the sum is across the 100 different dipole configurations
            % the sum is across the set of thresholds tested, so the number
            % returned reflects the hit rate based on the criterion that any of the
            % thresholds gave a hit.
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,1);bar(uvals{3}./100, Hx);
            ylim([0 1]); ylabel('hitrate coh 1-dip'); xlabel('coherence strength');
            legend(cellstr(num2str(phs)),'location','NOrthWest');
            dx = permute(cat(6, allCorig.D), [6 1 2 3 4 5]);
            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,2);bar(uvals{3}./100, Hx);
            ylim([0 1]); ylabel('hitrate coh 2-dip'); xlabel('coherence strength');
            dx = permute(cat(6, allCc.D), [6 1 2 3 4 5]);
            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,3);bar(uvals{3}./100, Hx);
            ylim([0 1]); ylabel('hitrate imag(coh)'); xlabel('coherence strength');
            dx = permute(cat(6, allC.D), [6 1 2 3 4 5]);
            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,4);bar(uvals{3}./100, Hx);
            ylim([0 1]); ylabel('hitrate coh 2-dip subsample'); xlabel('coherence strength');
            exportgraphics(gcf, sprintf('detection_rate_%d_%d_v2.pdf', uvals{1}(i1), uvals{2}(i2)));
            exportgraphics(gcf, sprintf('detection_rate_%d_%d_v2.png', uvals{1}(i1), uvals{2}(i2)));

        end
    end

    for i1 = 1:2
        for i2 = 1:3
            figure('position',[10 10 1120 840]);
            dx = permute(cat(6, allCorig1.D), [6 1 2 3 4 5]); % distance metric, goes to 0
            % when the simulated dipoles are within the 'blobs'

            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;

            % the sum is across the 100 different dipole configurations
            % the sum is across the set of thresholds tested, so the number
            % returned reflects the hit rate based on the criterion that any of the
            % thresholds gave a hit.
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,1);bar(uvals{3}./100, Hx-Hx(1,:));
            ylim([-.1 1]); ylabel('hitrate coh 1-dip'); xlabel('coherence strength');
            legend(cellstr(num2str(phs)),'location','NOrthWest');
            dx = permute(cat(6, allCorig.D), [6 1 2 3 4 5]);
            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,2);bar(uvals{3}./100, Hx-Hx(1,:));
            ylim([-.1 1]); ylabel('hitrate coh 2-dip'); xlabel('coherence strength');
            dx = permute(cat(6, allCc.D), [6 1 2 3 4 5]);
            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,3);bar(uvals{3}./100, Hx-Hx(1,:));
            ylim([-.1 1]); ylabel('hitrate imag(coh)'); xlabel('coherence strength');
            dx = permute(cat(6, allC.D), [6 1 2 3 4 5]);
            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            Hx = squeeze(sum(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02,6))>0, 1))./100;
            subplot(2,2,4);bar(uvals{3}./100, Hx-Hx(1,:));
            ylim([-.1 1]); ylabel('hitrate coh 2-dip subsample'); xlabel('coherence strength');
            exportgraphics(gcf, sprintf('detection_rate_%d_%d_v3.pdf', uvals{1}(i1), uvals{2}(i2)));
            exportgraphics(gcf, sprintf('detection_rate_%d_%d_v3.png', uvals{1}(i1), uvals{2}(i2)));

        end
    end
end

%% Plot False-positives

if(extra_figs)
    for i1 = 1:2  %#ok
        for i2 = 1:3
            figure('position',[10 10 1120 840]);
            dx = permute(cat(6, allCorig1.D), [6 1 2 3 4 5]);
            nx = permute(cat(6, allCorig1.N2), [6 1 2 3 4 5]);
            nx = squeeze(nx(:,i1,i2,:,:,3:end)-double(dx(:,i1,i2,:,:,3:end)<0.02));
            nx = reshape(nx, [100 40 7]);

            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            [Hx, ix] = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
            N = nan(size(ix)); for nn = 1:numel(ix), N(nn) = nanmean(nx(:,nn,ix(nn))); end

            subplot(2,2,1);bar(uvals{3}./100, N);
            ylim([0 max(N(:))]); ylabel('#FP coh 1-dip'); xlabel('coherence strength');
            legend(cellstr(num2str(phs)),'location','NOrthWest');

            dx = permute(cat(6, allCorig.D), [6 1 2 3 4 5]);
            nx = permute(cat(6, allCorig.N2), [6 1 2 3 4 5]);
            nx = squeeze(nx(:,i1,i2,:,:,3:end)-double(dx(:,i1,i2,:,:,3:end)<0.02));
            nx = reshape(nx, [100 40 7]);

            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            [Hx, ix] = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
            N = nan(size(ix)); for nn = 1:numel(ix), N(nn) = nanmean(nx(:,nn,ix(nn))); end

            subplot(2,2,2);bar(uvals{3}./100, N);
            ylim([0 max(N(:))]); ylabel('#FP coh 2-dip'); xlabel('coherence strength');

            dx = permute(cat(6, allCc.D), [6 1 2 3 4 5]);
            nx = permute(cat(6, allCc.N2), [6 1 2 3 4 5]);
            nx = squeeze(nx(:,i1,i2,:,:,3:end)-double(dx(:,i1,i2,:,:,3:end)<0.02));
            nx = reshape(nx, [100 40 7]);

            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            [Hx, ix] = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
            N = nan(size(ix)); for nn = 1:numel(ix), N(nn) = nanmean(nx(:,nn,ix(nn))); end

            subplot(2,2,3);bar(uvals{3}./100, N);
            ylim([0 max(N(:))]); ylabel('#FP imag(coh)'); xlabel('coherence strength');

            dx = permute(cat(6, allC.D), [6 1 2 3 4 5]);
            nx = permute(cat(6, allC.N2), [6 1 2 3 4 5]);
            nx = squeeze(nx(:,i1,i2,:,:,3:end)-double(dx(:,i1,i2,:,:,3:end)<0.02));
            nx = reshape(nx, [100 40 7]);

            %Hx = squeeze(sum(dx(:,2,2,:,:,7)<0.02))./100;
            [Hx, ix] = max(squeeze(sum(dx(:,i1,i2,:,:,3:end)<0.02))./100, [], 3);
            N = nan(size(ix)); for nn = 1:numel(ix), N(nn) = nanmean(nx(:,nn,ix(nn))); end

            subplot(2,2,4);bar(uvals{3}./100, N);
            ylim([0 max(N(:))]); ylabel('#FP 2-dip subsample'); xlabel('coherence strength');
            exportgraphics(gcf, sprintf('falsepositives_%d_%d.pdf', uvals{1}(i1), uvals{2}(i2)));
            exportgraphics(gcf, sprintf('falsepositives_%d_%d.png', uvals{1}(i1), uvals{2}(i2)));

        end
    end
end


%% Figure 7 and S4 - False positives


for i1 = snrs  % SNR
    for i2 = 1:3  % amplitudes
        figure('position',[10 10 840 600]); 
        hold all

        dx1 = permute(cat(6, allCorig1.D),  [6 1 2 3 4 5]);
        nx1 = permute(cat(6, allCorig1.N), [6 1 2 3 4 5]);
        nx1 = squeeze(nx1(:,i1,i2,2:end,:,1:end)-double(dx1(:,i1,i2,2:end,:,1:end)<0.02));
        hx1 = squeeze(double(dx1(:,i1,i2,2:end,:,1:end)<0.02));
        nnx1 = log10(squeeze(nanmedian(nanmedian(nanmedian(nx1)))));
        hhx1 = squeeze(mean(mean(mean(hx1))));

        dx2 = permute(cat(6, allCorig.D),  [6 1 2 3 4 5]);
        nx2 = permute(cat(6, allCorig.N), [6 1 2 3 4 5]);
        nx2 = squeeze(nx2(:,i1,i2,2:end,:,1:end)-double(dx2(:,i1,i2,2:end,:,1:end)<0.02));
        hx2 = squeeze(double(dx2(:,i1,i2,2:end,:,1:end)<0.02));
        nnx2 = log10(squeeze(nanmedian(nanmedian(nanmedian(nx2)))));
        hhx2 = squeeze(mean(mean(mean(hx2))));

        dx3 = permute(cat(6, allCc.D),  [6 1 2 3 4 5]);
        nx3 = permute(cat(6, allCc.N), [6 1 2 3 4 5]);
        nx3 = squeeze(nx3(:,i1,i2,2:end,:,1:end)-double(dx3(:,i1,i2,2:end,:,1:end)<0.02));
        hx3 = squeeze(double(dx3(:,i1,i2,2:end,:,1:end)<0.02));
        nnx3 = log10(squeeze(nanmedian(nanmedian(nanmedian(nx3)))));
        hhx3 = squeeze(mean(mean(mean(hx3))));

        dx4 = permute(cat(6, allC.D),  [6 1 2 3 4 5]);
        nx4 = permute(cat(6, allC.N), [6 1 2 3 4 5]);
        nx4 = squeeze(nx4(:,i1,i2,2:end,:,1:end)-double(dx4(:,i1,i2,2:end,:,1:end)<0.02));
        hx4 = squeeze(double(dx4(:,i1,i2,2:end,:,1:end)<0.02));
        nnx4 = log10(squeeze(nanmedian(nanmedian(nanmedian(nx4)))));
        hhx4 = squeeze(mean(mean(mean(hx4))));

        plot(nnx1, hhx1, 'd-', 'color', colors.single, 'linewidth', 2, ...
            'markerfacecolor', colors.single, 'markersize', 8);
        plot(nnx2, hhx2, 's-', 'color', colors.twodip, 'linewidth', 2, ...
            'markerfacecolor', colors.twodip, 'markersize', 8);
        plot(nnx3, hhx3, '+-', 'color', colors.imag, 'linewidth', 2, ...
            'markerfacecolor', colors.imag, 'markersize', 8);
        plot(nnx4, hhx4, 'o-', 'color', colors.subsam, 'linewidth', 2, ...
            'markerfacecolor', colors.subsam, 'markersize', 8);

        fontsize(gca, 14, 'pixels')
        xlabel('# of false positives', 'fontsize', 15);
        ylabel('hit rate', 'fontsize', 15);

        xlim([log10(2) log10(250)]);
        ylim([0 0.8]);
        set(gca, 'xtick', log10([2 5 10 20 50 100 200 250]));
        set(gca, 'xticklabel', [2 5 10 20 50 100 200 250]);
        set(gcf, 'color', 'w')

        for ii = 1:2  % no idea why, but this needs to be called twice to have an effect
            legend({'single dipole BF', 'two-dipole BF', 'single dipole BF, imag(coh)', 'two-dipole subsampling BF'}, ...
                'fontsize', 13, 'position', [455, 455, 340, 100], 'units', 'pixels', ...
                'autoupdate', 'off');
            legend boxoff
        end


        exportgraphics(gcf, sprintf('falsepositives_froc_%d_%d.pdf', uvals{1}(i1), uvals{2}(i2)));
        exportgraphics(gcf, sprintf('falsepositives_froc_%d_%d.png', uvals{1}(i1), uvals{2}(i2)), ...
            kwargs{:});

    end
end

%%

% % %%%%%%%%%%%%%%%%
% % load(fullfile('/project/3011020.09/jansch/tmp/results_noprewhiten_20210804_020', 'hits20210804_noprewhiten.mat'));
% % phs     = round(180*uvals{6}./(pi*100));
% % dipdist = squeeze(Dist(:,1,1,1,1,1,1,1,:));
% % sel     = dipdist(:,3)>15&dipdist(:,2)<5&dipdist(:,1)<5;
% %
% % figure;
% % Hx = max(squeeze(sum(Hits1(sel,:,2,:,:,:,:,:,3:end)))./sum(sel),[],3);
% % subplot(2,2,1);bar(uvals{5}./100, Hx);
% % ylim([0 1]); ylabel('hitrate coh 1-dip'); xlabel('coherence strength');
% % legend(cellstr(num2str(phs)),'location','NOrthWest');
% % Hx = max(squeeze(sum(HitsOrig(sel,:,2,:,:,:,:,:,3:end)))./sum(sel),[],3);
% % subplot(2,2,2);bar(uvals{5}./100, Hx);
% % ylim([0 1]); ylabel('hitrate coh 2-dip'); xlabel('coherence strength');
% % Hx = max(squeeze(sum(HitsC(sel,:,2,:,:,:,:,:,3:end)))./sum(sel),[],3);
% % subplot(2,2,3);bar(uvals{5}./100, Hx);
% % ylim([0 1]); ylabel('hitrate imag(coh)'); xlabel('coherence strength');
% % Hx = max(squeeze(sum(Hits(sel,:,2,:,:,:,:,:,3:end)))./sum(sel),[],3);
% % subplot(2,2,4);bar(uvals{5}./100, Hx);
% % ylim([0 1]); ylabel('hitrate coh 2-dip subsample'); xlabel('coherence strength');
% %
% % figure;
% % Hx = (squeeze(sum(Fp1(sel,:,2,:,:,:,:,:,8)))./sum(sel));
% % subplot(2,2,1);bar(uvals{5}./100, Hx);
% % ylim([0 15]); ylabel('#false positives coh 1-dip'); xlabel('coherence strength');
% % legend(cellstr(num2str(phs)),'location','NOrthWest');
% % Hx = (squeeze(sum(FpOrig(sel,:,2,:,:,:,:,:,8)))./sum(sel));
% % subplot(2,2,2);bar(uvals{5}./100, Hx);
% % ylim([0 15]); ylabel('#false positives coh 2-dip'); xlabel('coherence strength');
% % Hx = (squeeze(sum(FpC(sel,:,2,:,:,:,:,:,8)))./sum(sel));
% % subplot(2,2,3);bar(uvals{5}./100, Hx);
% % ylim([0 15]); ylabel('#false positives imag(coh)'); xlabel('coherence strength');
% % Hx = (squeeze(sum(Fp(sel,:,2,:,:,:,:,:,8)))./sum(sel));
% % subplot(2,2,4);bar(uvals{5}./100, Hx);
% % ylim([0 15]); ylabel('#false positives coh 2-dip subsample'); xlabel('coherence strength');
