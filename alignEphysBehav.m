function etho_spikes = alignEphysBehav(etho,sync_data,sp)

% extracts ethovision behavior output and checks that it aligns with ephys
% sync signal
% throughout script, e denotes ephys data, b denotes behav/etho data
% alj march 2018

%% get sync_data from neuropixels
sync = sync_data;
sync_on = sync.sync_on{1};
sync_trace = sync.sync_trace;
eTrials = length(sync.trial_on);
bTrials = length(etho.etho_photodiode);

if eTrials ~= bTrials
    disp('Photodiode & trials are mismatched.')
end

%% GET DATA FOR EACH TRIAL

for iTrial = 1:eTrials
    
    % sort out ephys data
    eTrialStart = sync.trial_on(iTrial);
    eTrialEnd = sync.trial_off(iTrial);
    eTrialLength = eTrialEnd - eTrialStart;

    ethoT = etho.recording_time{iTrial}; % get ethovision time
    
    % first we need to resample our ephys sync trace down to the ethovision
    % sampling rate
    trial_sync = sync.sync_trace(eTrialStart*2500:eTrialEnd*2500);
    nBins = length(etho.trial_time{iTrial}); % get time vector from the ethovision file
    % nBins = ceil(length(trial_sync)/62.5)
    binSize = floor(length(trial_sync)/nBins);
    for iBin = 1:nBins
        if iBin ~= nBins
            this_bin = trial_sync((binSize*(iBin-1)+1):binSize*iBin);
            this_bin_diff = diff([0;double(this_bin)]);
            if max(this_bin_diff) == 1 % if it turns on
                new_eSync(iBin) = 1;
            else
                new_eSync(iBin) = 0;
            end
        else
            this_bin = trial_sync((binSize*(iBin-1)+1):end);
            this_bin_diff = diff([0;double(this_bin)]);
            if max(this_bin_diff) == 1
                new_eSync(iBin) = 1;
            else
                new_eSync(iBin) = 0;
            end
        end
    end
    
    eSyncOn_IDs = find(new_eSync);
    eSyncOn_times = ethoT(eSyncOn_IDs);
    
    %     eSyncOn_IDs = find(sync_on >= eTrialStart & sync_on <= eTrialEnd); % get sync on times that happen during this period
    %     eSyncOn = sync_on(eSyncOn_IDs);
    % eSyncOn = eSyncOn - eSyncOn(1); % make it start at 0
    eSyncNum = length(eSyncOn_times);
    ePhotodiode = sync_data.photodiode(iTrial);
    
    % sort out behavioral data
    bTrialStart = etho.recording_time{iTrial}(1); % should always be zero, but making sure. Can also use Trialtime.
    bTrialEnd = etho.recording_time{iTrial}(end);
    bTrialLength = bTrialEnd - bTrialStart;
    bSyncTrace = etho.etho_sync{iTrial};
    bDiff = diff([0;double(bSyncTrace)]);
    bSyncOnIDs = find(bDiff<0); % it's reversed because etho vision is asking if the signal is low
    bSyncOn_times = ethoT(bSyncOnIDs);
    % bSyncOn = bSyncOn - bSyncOn(1); % make it start at 0 (otherwise it indexes at 1...?)
    bSyncNum = length(bSyncOn_times);
    
    % visually inspect sync traces
    figure;plot([eSyncOn_times eSyncOn_times],[0 1],'k')
    hold on; plot([bSyncOn_times bSyncOn_times],[0 1],'b')
    
    if eSyncNum ~= bSyncNum % if there are the same number of syncs, we're all good
        disp('Different sync traces, fix me?')
        disp(iTrial)
    end
    
    %% get relevant chunk of spikes
    
    %% get necessary info out of spike struct
    st = sp.st;
    cids = sp.cids;
    cgs = sp.cgs;
    clu = sp.clu;
    clusterDepths = sp.clusterDepths;
    
    %% count spikes in bins
    
    binSize = 1; %sec
    t = eTrialStart:binSize:eTrialEnd;
    numBins = length(t);
    
    inclCIDs = cids(cgs==2); % choose "good" clusters
    inclCIDdepths = clusterDepths(ismember(cids, inclCIDs)); % the depths of those clusters
    
    thisClu = clu(st>eTrialStart & st<eTrialEnd & ismember(clu,inclCIDs));
    thisST = st(st>eTrialStart & st<eTrialEnd & ismember(clu,inclCIDs));
    thisSTbins = ceil(thisST./binSize);
    
    % trick to do this binning with the sparse matrix creation function
    bin2d = full(sparse(double(thisClu+1), thisSTbins, 1));
    binnedSpikes = bin2d(inclCIDs+1,:); % then take just the rows according to the clusters you had
    [~, ii] = sort(inclCIDdepths);
    binnedSpikes = binnedSpikes(ii,:);
    
    binnedRates = binnedSpikes./binSize;
    w = gausswin(10);
    smoothedRates = filter(w,1,binnedRates);
    
    etho.binnedRates{iTrial} = binnedRates;
    etho.smoothedRates{iTrial} = smoothedRates;
    
end

etho_spikes = etho;


