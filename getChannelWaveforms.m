function theseWFs = getChannelWaveforms(filename,sp)

% function to extract channel waveforms
% you have to predefine the clusters & channels (i took these from
% kilosort)

% ALJ JULY 2018
% NP16 best channels = 109


%% file options
nTotalChan = 385;
chunkSize = 10000000;
nWFs = 200;
spikeWindow = 15;

%% which clusters/channels?
% clusters = [120 294 430 605 377 813 781 35 423 823 1349 792  ]; %NP16 (these are manually identified in phy)
% channels = [ 25  45 106  94 121 109 286 118 41 112  148 287  ]; %NP16
% clusters = [1205 1226 1183 1217 1261 1335 1312 1528 1242]; %NP14
% channels = [ 181  171  132  121  176  309   57   16  192]; %NP14
% clusters = [461  94 482 556 466 290 563 514 1011 525 513 176]; % not sure
% which animal this is
% channels = [375 361 312 302 297 229 213 208  193 181 181 172];

% or, we can pick channels based on certain characteristics
% here, we'll pick channels based on their depth & amplitude (how many
% spikes)
% this was used to pull out a putative common neuron throughout experiments in NP8

clusters = sp.cids(sp.clusterDepths>1500 & sp.clusterDepths < 2000 & sp.clusterAmps > 500)';
channels = [179:182];

nCh = length(channels);

%% open file

fid = fopen(filename, 'r');
dat = fread(fid, [nTotalChan chunkSize], '*int16');
fclose(fid);


%% get spike times for each channel & plot

for nClust = 1:length(clusters)
    
    figure;
    set(gcf,'Color','w','Position',[1         530        1291         275])
    
    for iCh = 1:nCh
        
        theseSpikeTimes = sp.ss(sp.clu == clusters(nClust));
 
        subplot(1,4,iCh);

        for i = 1:nWFs
            spikeStart = theseSpikeTimes(i) - spikeWindow;
            spikeEnd = theseSpikeTimes(i) + spikeWindow;
            theseWFs(i,:) = dat(channels(iCh),[spikeStart:spikeEnd]);
            plot(theseWFs(i,:),'Color',[.5 .5 .5])
            hold on
        end
                
        plot(mean(theseWFs),'LineWidth',2,'Color',[1 .8 0]);
        title(num2str(channels(iCh)));
        axis tight
        box off
        allWFs{iCh} = theseWFs;
        clear theseWFs
    end  
end





