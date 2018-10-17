function massiveNeuropixScript

% largely based on nicks code, cortex http://data.cortexlab.net/dualPhase3/data/script_dualPhase3.m

%% Script to do some basic loading of data and plotting

% make sure npy-matlab and spikes are in matlab path
addpath(genpath('C:\Users\Ashley\Documents\npy-matlab'))
addpath(genpath('C:\Users\Ashley\Documents\spikes'))

% add data dir
if ismac
    datadir = '/Users/ajuavinett/Desktop/Neuropixels/12062017/';
else
    datadir = 'C:\Neuropixels\NP6\12062017';
end

addpath(genpath(datadir))

%% load some spikes and compute some basic things 
 
% This contains the x- and y-coordinates in space of every recording site
% as "xcoords" and "ycoords" (this information is also in
% channel_positions.npy) and the list of which channels are connected (only
% 374 out of the 384 can record data; "connected" is logical, size 384x1)
if ismac
    load('/Volumes/churchland_hpc_home/ajuavine/Neuropixels/neuropixPhase3_kilosortChanMap_opt4.mat')
else
    load(fullfile('C:\Users\Ashley\Documents\neuropixPhase3_kilosortChanMap_opt4.mat'))
end

yc = ycoords(connected); xc = xcoords(connected);
Fs = 30000;

% clu is a length nSpikes vector with the cluster identities of every spike
clu = readNPY(fullfile(datadir,'spike_clusters.npy'));

% ss is a length nSpikes vector with the spike time of every spike (in samples)
ss = readNPY(fullfile(datadir, 'spike_times.npy'));

% convert to times in seconds
st = double(ss)/Fs;

% spikeTemplates is like clu, except with the template numbers rather than cluster numbers.
% Each spike was extracted by one particular template(identified here)
% but when templates were merged in the manual sorting, the spikes of both take on a new cluster identity in clu.
% So spikeTemplates reflects the original output of the algorithm; clu is the result of manual sorting. 
spikeTemplates = readNPY(fullfile(datadir,'spike_templates.npy')); % note: zero-indexed

% tempScalingAmps is a length nSpikes vector with the "amplitudes":
% each spike is extracted by a particular template after scaling the
% template by some factor - that factor is the amplitude. The actual
% amplitude of the spike is this amplitude multiplied by the size of the
% template itself - we compute these later. 
tempScalingAmps = readNPY(fullfile(datadir, 'amplitudes.npy'));

% get cluster groups (see discussion in http://data.cortexlab.net/dualPhase3/data/script_dualPhase3.m)
[cids, cgs] = readClusterGroupsCSV(fullfile(datadir,'cluster_groups.csv'));

% find and discard spikes corresponding to noise clusters
noiseClusters = cids(cgs==0);

st = st(~ismember(clu, noiseClusters));
ss = ss(~ismember(clu, noiseClusters));
spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
clu = clu(~ismember(clu, noiseClusters));
cgs = cgs(~ismember(cids, noiseClusters));
cids = cids(~ismember(cids, noiseClusters));

% temps are the actual template waveforms. It is nTemplates x nTimePoints x
% nChannels (in this case 1536 x 82 x 374). These should be basically
% identical to the mean waveforms of each template
temps = readNPY(fullfile(datadir,'templates.npy'));

% The templates are whitened; we will use this to unwhiten them into raw
% data space for more accurate measurement of spike amplitudes; you would
% also want to do the same for spike widths. 
winv = readNPY(fullfile(datadir,  'whitening_mat_inv.npy'));

% compute some more things about spikes and templates; see function for
% documentation
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW] = ...
    templatePositionsAmplitudes(temps, winv, yc, spikeTemplates, tempScalingAmps);

% convert to uV according to the gain and properties of the probe
% 0.6 is the range of voltages acquired (-0.6 to +0.6)
% 512 is the bit range (-512 to +512, 10bits)
% 500 is the gain factor I recorded with 
% 1e6 converts from volts to uV
spikeAmps = spikeAmps*0.6/512/500*1e6; 

% save everything in a struct
sp.name = datadir;
sp.clu = clu;
sp.ss = ss;
sp.st = st;
sp.spikeTemplates = spikeTemplates;
sp.tempScalingAmps = tempScalingAmps;
sp.cgs = cgs;
sp.cids = cids;
sp.yc = yc;
sp.xc = xc;
sp.ycoords = ycoords;
sp.xcoords = xcoords;
sp.temps = temps;
sp.spikeAmps = spikeAmps;
sp.templateYpos = templateYpos;
sp.tempAmps = tempAmps;
sp.spikeDepths = spikeDepths;
sp.tempsUnW = tempsUnW;


%% depths and amplitudes of clusters (as the mean depth and amplitude of all of their constituent spikes)
% get firing rates here also
recordingDuration = sp.st(end)-sp.st(1);
sd = sp.spikeDepths;
sa = sp.spikeAmps;
% using a super-tricky algorithm for this - when you make a sparse
% array, the values of any duplicate indices are added. So this is the
% fastest way I know to make the sum of the entries of sd for each of
% the unique entries of clu
[cids, spikeCounts] = countUnique(clu);
q = full(sparse(double(clu+1), ones(size(clu)), sd));
q = q(cids+1);
clusterDepths = q./spikeCounts; % had sums, so dividing by spike counts gives the mean depth of each cluster

q = full(sparse(double(clu+1), ones(size(clu)), sa));
q = q(cids+1);
clusterAmps = q./spikeCounts;

sp.clusterDepths = clusterDepths';
sp.clusterAmps = clusterAmps';
sp.firingRates = spikeCounts'./recordingDuration;


%% basic plot of clusters over depth 
% on each probe, higher depth numbers are superficial, i.e. nearer to the
% top of the brain; lower numbers are deeper, nearer the tip

v1Borders = [2797 3840]; % determined by manual inspection
hcBorders = [1634 2797];
thalBorders = [0 1634];

f = figure;

cd = sp.clusterDepths;
ca = sp.clusterAmps;
cgs = sp.cgs;

xx = rand(size(cgs));

scatter(xx,cd,ca/5);
hold on;
scatter(xx(cgs==2),cd(cgs==2),ca(cgs==2)/5)
title([sp.name ' probe, neuron depths and amplitudes']);
xlabel('random value for visualization')
ylabel('depth on probe (um)')
legend({'MUA/Unsorted', 'Good'});
plot([0 1], v1Borders(1)*[1 1], 'k--', 'LineWidth', 2.0);
plot([0 1], hcBorders(1)*[1 1], 'k--', 'LineWidth', 2.0);

v1Count = sum(sp.clusterDepths>=v1Borders(1) & sp.cgs==2);
hcCount = sum(sp.clusterDepths>=hcBorders(1) & sp.clusterDepths<hcBorders(2) & sp.cgs==2);
thalCount = sum(sp.clusterDepths<thalBorders(2) & sp.cgs==2);


% fprintf(1, 'recorded %d v1 neurons, %d hippocampus, %d thalamus\n', v1Count, hcCount, thalCount);

%% raster plot

win = [1306 1308];
rasterScale = 50;

f = figure; set(f, 'Color', 'w');

st = sp.st;
inclSpikes = st>win(1) & st<=win(2);
st = st(inclSpikes);
clu = sp.clu; clu = clu(inclSpikes);
cids = sp.cids;
cgs = sp.cgs;
sd = sp.spikeDepths; sd = sd(inclSpikes);
    
co = get(gca, 'ColorOrder'); nc = size(co,1);

for c = 1:length(sp.cids)
    if sp.cgs(c)==2
        thisColor = co(mod(c,nc)+1,:);
    else
        thisColor = [0.7 0.7 0.7]; % grey for unsorted and mua
    end
    
    [xx, yy] = rasterize(st(clu==cids(c)));
    theseSD = sd(clu==cids(c));
    yy(1:3:end) = yy(1:3:end)+theseSD';
    yy(2:3:end) = yy(2:3:end)*rasterScale+theseSD';
    
    plot(xx,yy, 'Color', thisColor, 'LineWidth', 1.5);
    hold on;
end
    
box off;
% set(gca, 'XTick', []);
%ylabel('depth on probe (um)')
%xlabel('time (sec)');

%% plot some stimulus responses
sync_data = masterSyncChan(datadir);
stimStarts = sync.photodiode;
stimStarts = stimStarts(find(stimStarts > sync.center(1)));
% stimIDs = ones(length(stimStarts),1);
% stimIDs = [1;2;1;1;2;2;3;3;3;1;1;2;2;3;3;3;1;2;2;3;1]; %for NP6_run3
% stimIDs = [2;1;3]; %for NP7
% stimIDs = [2 4 3 1 3 2 1 4 2 4 3 1 1 3 2]; % NP6_005
% stimIDs = [3 3 1 1 2 1 2 2 1 3 3 1 1 3 3 1 2 1 2 1 1 3 2 2 1 1 1 2 2 1 1 3 3 1 1 2 3 3 1 2]; %NP6_headfixed2

st = sp.st;
clu = sp.clu;
cgs = sp.cgs;
cids = sp.cids;
clusterDepths = sp.clusterDepths;
FRs = sp.firingRates;

v1Borders = [2797 3840]; % determined by manual inspection

% select Good neurons in V1 with at least some reasonable spike rate
%inclClusters = cids(cgs==2 & clusterDepths>v1Borders(1) & clusterDepths<=v1Borders(2) & FRs>0.5); 
inclClusters = cids(cgs==2 & FRs>1); 

window = [-0.5 2.5]; % time window over which to make a psth

% use arrow keys to move through clusters, see function help for other
% controls. E.g. try pressing "t". 
%
% Flip through at least until you see cluster 657, which is a good example
psthViewer(st(ismember(clu,inclClusters)), clu(ismember(clu, inclClusters)), stimStarts, window, stimIDs);

%% count spikes in bins

binSize = 1; %sec

t = 0:binSize:st(end); % whole recording
nBins = length(t);

inclCIDs = cids(cgs==2); % choose "good" clusters
inclCIDdepths = clusterDepths(ismember(cids, inclCIDs)); % the depths of those clusters

thisClu = clu(st>0 & ismember(clu,inclCIDs));
thisST = st(st>0 & ismember(clu,inclCIDs));
thisSTbins = ceil(thisST./binSize);

% trick to do this binning with the sparse matrix creation function
bin2d = full(sparse(double(thisClu+1), thisSTbins, 1));
binnedSpikes = bin2d(inclCIDs+1,:); % then take just the rows according to the clusters you had
[sortedCIDdepths, ii] = sort(inclCIDdepths);
binnedSpikes = binnedSpikes(ii,:);

binnedRates = binnedSpikes./binSize;

%% simple visualization, with colormap as spike rate, x-axis time, y-axis cluster index

tstart = 0; 
winSize = 60; % seconds

figure; 
nSamps = sum(t>tstart & t<tstart+winSize);
im = imagesc((1:nSamps)*binSize, 1:size(binnedRates,1), binnedRates(:, t>tstart & t<tstart+winSize));
set(gca, 'YDir', 'normal');
caxis([0 70])
colorbar
xlabel('time (sec)')
ylabel('<-- thalamus, hc, vis ctx -->');
title(sprintf('start time %d seconds', tstart))

%% play through the recording as a movie

while 1
    set(im, 'CData', binnedRates(:, t>tstart & t<tstart+winSize));     
    title(sprintf('start time %d seconds', tstart))
    drawnow; 
    tstart = tstart+1; 
end

%% correlation between the binned spike rates

corrMat = corr(binnedRates');

%% plot correlation

v1Borders = [2797 3840]; % determined by manual inspection
hcBorders = [1634 2797];
thalBorders = [0 1634];

figure; 
imagesc(corrMat-diag(diag(corrMat)));
colormap(colormap_blueblackred)
axis square
caxis([-0.6 0.6])
set(gca, 'YDir', 'normal');
firstHC = find(sortedCIDdepths>hcBorders(1),1);
firstV1 = find(sortedCIDdepths>v1Borders(1),1);
hold on;
lims = ylim();
plot(lims, [firstHC firstHC], 'w--');
plot([firstHC firstHC], lims, 'w--');
plot(lims, [firstV1 firstV1], 'w--');
plot([firstV1 firstV1], lims, 'w--');

