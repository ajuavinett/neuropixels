function plotStimResponses(sp,sync)

% to plot raw responses to stimuli

%% plot some stimulus responses
stimStarts = sync.photodiode;
% stimStarts = stimStarts(8:47); % this is necessary bc there are pulses before actual start in /Users/ajuavinett/Desktop/Neuropixels/12062017/
% stimStarts = stimStarts(find(stimStarts > sync.center(1)));
% stimIDs = ones(length(stimStarts),1);
% stimIDs = [1;2;1;1;2;2;3;3;3;1;1;2;2;3;3;3;1;2;2;3;1]; %for NP6_run3
% stimIDs = [2;1;3]; %for NP7
% stimIDs = [2 4 3 1 3 2 1 4 2 4 3 1 1 3 2]; % NP6_005
% stimIDs = [3 3 1 1 2 1 2 2 1 3 3 1 1 3 3 1 2 1 2 1 1 3 2 2 1 1 1 2 2 1 1 3 3 1 1 2 3 3 1 2]; %NP6_headfixed2

stimIDs = makeUniqueStimID(sync.contrast,sync.sound_bit,sync.type);

st = sp.st;
clu = sp.clu;
cgs = sp.cgs;
cids = sp.cids;
clusterDepths = sp.clusterDepths;
FRs = sp.firingRates;

v1Borders = [2797 3840]; % determined by manual inspection

%% select neurons based on some criteria
%inclClusters = cids(cgs==2 & clusterDepths>v1Borders(1) & clusterDepths<=v1Borders(2) & FRs>0.5); 
% inclClusters = cids(cgs==2 & FRs>1); 
inclClusters = cids(cgs==2);

window = [-0.5 2.5]; % time window over which to make a psth

psthViewer(st(ismember(clu,inclClusters)), clu(ismember(clu, inclClusters)), stimStarts, window, stimIDs);
