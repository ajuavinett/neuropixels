function [psthSm, stderr] =  generatePSTHforGroup(baSm,bins,myData,group,grouplabels)

% baSm is smoothed trial x time matrix for neuron
% group is your grouping variable

nGroups = length(grouplabels);
psthSm = zeros(nGroups, numel(bins));

if myData.params.showErrorShading
    stderr = zeros(nGroups, numel(bins));
end
for g = 1:nGroups
    this_hist = baSm(group==grouplabels(g),:);
    [nRepeats,~] = size(this_hist); % added this for trials with only 1 repeat
    if nRepeats ~= 1
        psthSm(g,:) = mean(this_hist);
    else
        psthSm(g,:) = this_hist;
    end
    if myData.params.showErrorShading && nRepeats > 1
        stderr(g,:) = std(baSm)./sqrt(size(baSm,1));
    end
end