function masterNeuropix(expt_name)

%% READ ME
% This is the master script to analyze neuropixels experiments. It creates
% spikes.mat, sync.mat, etho.mat, and meta.mat structures, plots depth of
% neurons and a sample raster, and creates panels for each good cluster.
% You can also run this script to load those structs once they have been made.

% INPUTS: expt_name = 'NP8_001'; this works for both freely moving and
% headfixed.

% IN ORDER TO RUN THIS, you need to have the outputs of kilosort (not
% including the huge .bin file) in a folder. The script "getDataDir" should
% point to this folder.
% You should also move the lf.bin file to the same folder - this will be
% used to extract the sync pulses from the Neuropixels acquisition.

% FOR FREELY MOVING EXPERIMENTS
%   * move the trialArray for this experiment into the same folder. It will be called in "makeStimFile." 
%   * make sure Ethovision Raw output is already exported. It will be called in "readEthovisionRaw.

% NOTE ABOUT SYNC PULSES
%   * very often, there are stray photodiode pulses. at the moment, these
%   are manually removed by plotting the sync pulse (plotSyncTrace) and
%   looking for them.

%% SCRIPT INPUTS

if nargin == 0
    datadir = uigetdir('Pick a data dir');
    fileparts = strsplit(datadir,'\'); 
    expt_name = fileparts{3};
else
    datadir = getDataDir(expt_name);
end

addpath(genpath(datadir)) % make sure datadir is a matlab path
[meta, sp] = makeMetaSpikes(datadir,expt_name); % make/load meta & spikes files

%% GET SYNC FILE (OR MAKE IT)
sync_file = fullfile(datadir,strcat(expt_name,'_','sync.mat'));

if ~exist(sync_file,'file') % if we didn't already save a sync file
    disp('making a sync_data file.')
    sync_data = masterSyncChan(1); % make sync file
    switch meta.expt_type
        case 'free'
            plotSyncTraces(sync_data,1)% visually inspect sync file    
            sync_data = makeStimFile(datadir,sync_data);
            save(sync_file,'sync_data');
        case 'headfixed'
            % plotSyncTraces(sync_data,0)% visually inspect sync file                
            edited_sync_data = makeHeadfixedSync(expt_name,sync_data,meta.stim_type);
            sync_data.photodiode = edited_sync_data.photodiode;
            sync_data.stimIDs = edited_sync_data.stimIDs;
            sync_data.trialPulses = ones(length(sync_data.stimIDs),1);
            save(sync_file,'sync_data')
    end
else
    load(sync_file)
end

%% BRING IN ETHOVISION FILE
if strcmp(meta.expt_type,'free')
    etho_file = strcat(datadir,'\',expt_name,'_etho.mat');
    if ~exist(etho_file,'file')
        disp('creating ethovision .mat file')
        etho = readEthovisionRaw(expt_name,sync_data);
    else
        disp('loading ethovision file')
        load(etho_file);
    end
    % etho_spikes = alignEphysBehav(etho,sync_data,sp);
    % plotBehaviorWithSpikes(etho_spikes);
end

%% ANALYZE MOVEMENT
switch meta.expt_type
    case 'headfixed'
        mvmtNeurons(sp,sync_data);
end

%% PLOT EXPERIMENT FIGURES
overviewFig = figure; set(overviewFig, 'Color', 'w','Position',[1   574   609   420]);
subplot(121); plotClustersbyDepth(meta,sp);
subplot(122); plotNeuropixRaster(meta,sp,[150 200],sync_data); % meta, sp, time window, save_bit

savedir = fullfile(meta.datadir,'Figures');
if ~exist(savedir,'dir'); mkdir(savedir);end
saveas(gcf,strcat(savedir,'\','neuronsbydepth'),'fig')
saveas(gcf,strcat(savedir,'\','neuronsbydepth'),'svg')

%% PLOT FIGURE BY CLUSTER
switch meta.expt_type
    case 'free'
        clusterViewer(sp,sync_data,meta,etho)
    case 'headfixed'
        switch meta.stim_type
            case 'AVLR'
                stimViewer(sp,sync_data,meta)
                % plotStimResponses(sp,sync_data); % old
            case 'ori'
                stimViewer(sp,sync_data,meta)
                % plotOriResponses(sp,sync_data); % old
            case 'retinotopy'
                plotRetResponses(sp,sync_data);
        end
end

