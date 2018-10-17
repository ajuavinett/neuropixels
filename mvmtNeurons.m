function [binnedRates,velocity,inclCIDs] = mvmtNeurons(sp,sync,plot_bit)

% script to get rough mvmt trace on wheel
% and compare this to neural activity of isolated neurons

mvmt1_trace = sync.mvmt1_trace; % get first mvmt trace
mvmt2_trace = sync.mvmt2_trace;
mvmt1_change = abs(diff(mvmt1_trace)); % did it change?
IDs = find(mvmt1_change);

mvmt = zeros(length(mvmt1_trace),1);

for i = 1:length(IDs)
    thisID = IDs(i);
    if mvmt1_trace(thisID) == mvmt2_trace(thisID)
        mvmt(thisID) = 1; % originally these were flipped but it seemed like it was all backward?
    else
        mvmt(thisID) = -1;
    end
end

mvmt_sum = cumsum(mvmt);
mvmt_sum_sec = downsample(mvmt_sum,2500);

samples = 5000; % number of samples to average for velocity
time = (samples/2500)*1000;
sample_vector = [1:samples:length(mvmt_sum)];

change = [];
velocity = [];
for i = 1:(length(sample_vector)-1)
    change(i) = mvmt_sum(sample_vector(i+1)) - mvmt_sum(sample_vector(i));
    velocity(i) = change(i)/time;
end

% old way of getting crude velocity
expt_length = mvmt(end);
mvmt_bins = ceil(length(mvmt)); % for 1 ms bins
[mvmt_hist,edges] = hist(mvmt,mvmt_bins);

%% get necessary info out of struct
st = sp.st;
cids = sp.cids;
cgs = sp.cgs;
clu = sp.clu;
clusterDepths = sp.clusterDepths;

%% count spikes in bins

binSize = time/1000; %sec

t = 0:binSize:st(end); % whole recording
nBins = length(t)-1;

inclCIDs = cids; % take all
% inclCIDs = cids(cgs==2); % choose "good" clusters
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

%% plot!
if plot_bit
    plotMvmtNeurons(binnedRates,velocity,inclCIDs)
end