function plotOriResponses(sp,sync)

% to plot raw responses to drifting gratings stimuli

%% plot some stimulus responses
stimStarts = sync.photodiode;
stimIDs = sync.stimIDs;

[nTrials,nVars] = size(stimIDs);

if nVars > 1
        % find all possible combinations
        for i = 1:nTrials
            new_stimIDs{i} = strcat(num2str(stimIDs(i,1)),num2str(stimIDs(i,2)));
        end
        
        stimIDs = new_stimIDs;
end
    
st = sp.st;
clu = sp.clu;
cgs = sp.cgs;
cids = sp.cids;
clusterDepths = sp.clusterDepths;
FRs = sp.firingRates;

% v1Borders = [2797 3840]; % determined by manual inspection

%% select neurons based on some criteria
%inclClusters = cids(cgs==2 & clusterDepths>v1Borders(1) & clusterDepths<=v1Borders(2) & FRs>0.5); 
inclClusters = cids(cgs==2 & FRs>1); 

window = [-0.5 2.5]; % time window over which to make a psth

retViewer(st(ismember(clu,inclClusters)), clu(ismember(clu, inclClusters)), stimStarts, window, stimIDs);
