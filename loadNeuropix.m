function [meta, sp] = loadNeuropix(datadir,expt_name)

% script to load neuropixels data into matlab for various purposes

%% INPUTS, OUTPUTS & NECESSARY DIRECTORIES
% datadir is the directory that contains the .bin files as well as kilosort output
% needs spikes and npy-matlab directories
% outputs sp struct which has what you need
% originally by n. steinmetz, edited by a. juavinett

%% SET UP SCRIPT
addpath(genpath(datadir)) % make sure datadir is a matlab path
meta = {}; % set up our empty structs
sp = {};

%% META DATA
Fs = 30000; % sampling rate
[expt_parts, ~] = strsplit(expt_name,'_');
anim = expt_parts(1);
anim = anim{1};

%% GET OPTION, CHANNEL MAP, BORDERS
% channel_positions.npy contains the x- and y-coordinates in space of every recording site
% we can also do this by directly loading channel_positions.py, but xcoords
% and ycoords are only connected ch.

if strcmp(anim,'NP6')
    borders = [3000 3840; 0 3000]; %all midbrain, all day.
    probe_option = 4;
elseif strcmp(anim,'NP8')
    borders = [3480 3840; 3280 3480; 2380 3280; 0 2380]; %v1, commissure, subiculum, mb
    probe_option = 1;
elseif strcmp(anim,'NP9')
    borders = [3000 3840; 0 3000]; %all midbrain, all day.
    probe_option = 1;
elseif strcmp(anim,'NP16')
    borders = [4780        6260;   3980        4780;      3200        3980;        2420        3200];
    probe_option = 3;
elseif strcmp(anim,'NP14') || strcmp(anim,'NP15')
    borders = [3000 3840; 0 3000]; %need to adjust this once I look... (update: i manually fixed these after creating meta struct)
    probe_option = 3;
end

if isstring(borders) % because the inputdlg is a string...
    borders = str2num(borders{1});
end

if isstring(probe_option)
    probe_option = str2num(probe_option{1});
end

%% LOAD THE CORRECT CHANNEL MAP
if strcmp(anim,'NP8') || strcmp(anim,'NP9') || strcmp(expt_name,'NP14_000') || strcmp(expt_name,'NP14_001')
    [~, ~, xcoords, ycoords, connected, ~] = makeForPRBimecP3(probe_option);
elseif strcmp(expt_name,'NP14_003') || strcmp(expt_name,'NP14_004') || strcmp(expt_name,'NP14_005') || strcmp(expt_name,'NP14_006') || strcmp(expt_name,'NP14_007') || strcmp(expt_name,'NP14_008') || strcmp(expt_name,'NP14_009') || strcmp(expt_name,'NP14_010') || strcmp(expt_name,'NP14_011')
    load('A:\data\ashley_looming\Neuropixels\ChannelMaps\neuropixPhase3_kilosortChanMap_opt3_0-100_bank1.mat')
elseif strcmp(anim,'NP15')
    load('A:\data\ashley_looming\Neuropixels\ChannelMaps\neuropixPhase3_kilosortChanMap_opt3_0-140_bank1.mat')
elseif strcmp(anim,'NP16')
    load('A:\data\ashley_looming\Neuropixels\ChannelMaps\neuropixPhase3_kilosortChanMap_opt3_0-240_bank1.mat')
end

    xc = xcoords(connected); yc = ycoords(connected); % create different vectors for only connected ch


%% SAVE META DATA
meta.datadir = datadir;
meta.name = expt_name;
meta.anim = anim;
meta_filename = fullfile(datadir,strcat(expt_name,'_','meta.mat'));
meta.timestamp = datetime;
[meta.expt_type,meta.stim_type] = getExptStimType(datadir);
meta.borders = borders;
meta.probe_option = probe_option;
save(meta_filename,'meta');

%% LOAD SPIKES & COMPUTE SOME THINGS

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
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(temps, winv, yc, spikeTemplates, tempScalingAmps);

% convert to uV according to the gain and properties of the probe
% 0.6 is the range of voltages acquired (-0.6 to +0.6)
% 512 is the bit range (-512 to +512, 10bits)
% 500 is the gain factor I recorded with
% 1e6 converts from volts to uV
spikeAmps = spikeAmps*0.6/512/500*1e6;

%% ONLY KEEP THE GOOD/MUA SHIT
goodClusters = cids(cgs==2);
muaClusters = cids(cgs==1);
myClusters = [goodClusters,muaClusters]; % we are taking both for now

% below are all the variables that are spikes length
st = st(ismember(clu, myClusters));
ss = ss(ismember(clu, myClusters));
spikeTemplates = spikeTemplates(ismember(clu, myClusters));
tempScalingAmps = tempScalingAmps(ismember(clu, myClusters));
spikeAmps = spikeAmps(ismember(clu, myClusters));
spikeDepths = spikeDepths(ismember(clu, myClusters));
clu = clu(ismember(clu, myClusters));

% these are the variables that are clusters length
cgs = cgs(ismember(cids, myClusters));
cids = cids(ismember(cids, myClusters));

% currently not dealing with template variables

%% DEPTHS & AMPLITUDES OF CLUSTERS
% as the mean depth and amplitude of all of their constituent spikes)
% get firing rates here also
recordingDuration = st(end)-st(1);
% using a super-tricky algorithm for this - when you make a sparse
% array, the values of any duplicate indices are added. So this is the
% fastest way I know to make the sum of the entries of sd for each of
% the unique entries of clu
[cids, spikeCounts] = countUnique(clu);
q = full(sparse(double(clu+1), ones(size(clu)),spikeDepths));
q = q(cids+1);
clusterDepths = q./spikeCounts; % had sums, so dividing by spike counts gives the mean depth of each cluster

q = full(sparse(double(clu+1), ones(size(clu)),spikeAmps));
q = q(cids+1);
clusterAmps = q./spikeCounts;

clear q

%% SAVE SP STRUCT
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
sp.waveforms = waveforms;
sp.spikeDepths = spikeDepths;
sp.tempsUnW = tempsUnW;
sp.clusterDepths = clusterDepths';
sp.clusterAmps = clusterAmps';
sp.firingRates = spikeCounts'./recordingDuration;

spikes_file = fullfile(datadir,strcat(expt_name,'_','spikes.mat'));
save(spikes_file,'sp');
