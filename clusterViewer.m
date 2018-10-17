function clusterViewer(sp,sync_data,meta,etho)
%
% Controls:
% - c: dialog box to pick a new cluster ID number
% - t: toggle showing psth traces for each grouping variable or just the
% overall. If showing just overall, raster is sorted chronologically. If
% showing by grouping variable, raster is sorted by that variable.
% - r: select a new range within which to count spikes for the tuning curve

%% SET UP FIGURE
f = figure;
set(gcf,'Position',[1  39  1366 634],'Color','w')
expt_panel = uipanel('Title','Expt Info','FontSize',12,'BackgroundColor','white',...
    'Position',[.05 .65 .2 .3]);
id_panel = uipanel('Title','Cluster Info','FontSize',12,'BackgroundColor','white',...
    'Position',[.05 .25 .2 .3]);
wf_panel= uipanel('Parent',id_panel,'Title','Template','FontSize',12,...
    'BackgroundColor','w','Position',[.5 .1 .4 .9]);
wf_ax = axes('parent',wf_panel);
psth_panel = uipanel('Title','PSTH','FontSize',12,'BackgroundColor','white',...
    'Position',[.3 .65 .3 .3]);
psth_ax = axes('parent',psth_panel);
raster_panel = uipanel('Title','Raster Plot','FontSize',12,'BackgroundColor','white',...
    'Position',[.3 .35 .3 .3]);
raster_ax = axes('parent',raster_panel);
tuning_panel = uipanel('Title','Tuning','FontSize',12,'BackgroundColor','w',...
    'Position',[.3 .05 .3 .3]);
tuning_ax = axes('parent',tuning_panel);
trialFR_panel = uipanel('Title','Trial Firing Rate','FontSize',12,'BackgroundColor','white',...
    'Position',[.65 .05 .3 .3]);
trialFR_ax = axes('parent',trialFR_panel);
velocity_panel = uipanel('Title','Velocity','FontSize',12,'BackgroundColor','white',...
    'Position',[.65 .35 .3 .3]);
velocity_ax = axes('parent',velocity_panel);
behav_panel = uipanel('Title','Behavior','FontSize',12,'BackgroundColor','white',...
    'Position',[.65 .65 .3 .3]);
behav_ax = axes('parent',behav_panel);

%stability_panel = uipanel('Title','Stability','FontSize',12,'BackgroundColor','white','Position',[.05 .05 .9 .2]);

%% SET PARAMS
params.smoothSize = 20; % in msec, stdev of gaussian smoothing filter
params.clusterIndex = 1;
params.trialIndex = 1;
params.rasterScale = 3; % height of ticks
params.window = [-0.5 2]; % time window over which to make a stimulus psth
params.startRange = params.window(1);
params.stopRange = params.window(2);
params.showAllTraces = true;
params.showErrorShading = true;
params.binSize = 0.001;

%% ORGANIZE MY DATA (myData is what gets updated when you click through the figure)
myData.spikeTimes = sp.st;
myData.clu = sp.clu;
myData.spikeDepths = sp.spikeDepths;
myData.clusterDepths = sp.clusterDepths;
myData.cids = sp.cids;
myData.cgs = sp.cgs;
myData.spikeTemplates = sp.spikeTemplates;
myData.waveforms = sp.waveforms;
myData.borders = meta.borders;
myData.sync_data = sync_data;
myData.stim_type = meta.stim_type;
myData.expt_type = meta.expt_type;

eventTimes = sync_data.photodiode;
if ~issorted(eventTimes)
    [eventTimes, ii] = sort(eventTimes);
    trGroups = trGroups(ii);
end

%% THIS IS WHERE YOU DECIDE WHAT GETS GROUPED, AND WHICH STIM TIMES YOU USE
if median(sync_data.trialPulses) == 3
    myData.eventTimes = eventTimes(1:3:end); % take every third event time
    if isfield(sync_data,'stimIDs')
        myData.trGroups = sync_data.stimIDs(1:3:end);
    else
        myData.sync_data.stimIDs = makeUniqueStimID(sync_data.contrast,sync_data.sound_bit,sync_data.type);
        myData.trGroups = myData.sync_data.stimIDs(1:3:end);
    end
else
    myData.eventTimes = eventTimes(:);
    if isfield(sync_data,'stimIDs')
        myData.trGroups = sync_data.stimIDs;
    else
        myData.sync_data.stimIDs = makeUniqueStimID(sync_data.contrast,sync_data.sound_bit,sync_data.type);
        myData.trGroups = myData.sync_data.stimIDs;
    end
end

myData.clusterIDs = unique(myData.clu);
myData.trGroupLabels = unique(myData.trGroups);
myData.nGroups = length(myData.trGroupLabels);
myData.plotAxes = [psth_ax,raster_ax,tuning_ax,wf_ax,behav_ax,velocity_ax,trialFR_ax];
myData.panels = [id_panel,wf_panel,behav_panel,expt_panel,psth_panel];
myData.params = params;
if nargin > 3; myData.etho = etho; end
myData.meta = meta;

set(f, 'UserData', myData);
set(f, 'KeyPressFcn', @(f,k)psthViewerCallback(f, k));

psthViewerPlot(f)
end

function psthViewerPlot(f)

myData = get(f,'UserData');

%% ORGANIZE STIMULUS DRIVEN DATA AND COMPUTE CURVES
% pick the right spikes
st = myData.spikeTimes(myData.clu==myData.clusterIDs(myData.params.clusterIndex));

% compute everything
%[psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st, myData.eventTimes, myData.params.window, 0.001);
[psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(st,myData.eventTimes, myData.params.window, myData.params.binSize);
trGroupLabels = myData.trGroupLabels;
nGroups = myData.nGroups;
inclRange = bins>0 & bins<=1; % for spike counts (for tuning curve) only, just take one second
spikeCounts = sum(ba(:,inclRange),2)./(myData.params.stopRange-myData.params.startRange);

% PSTH smoothing filter
gw = gausswin(round(myData.params.smoothSize*6),3);
smWin = gw./sum(gw);

% smooth ba
baSm = conv2(smWin,1,ba', 'same')'./myData.params.binSize;

% construct psth(s) and smooth it
if myData.params.showAllTraces
    psthSm = zeros(nGroups, numel(bins));
    if myData.params.showErrorShading
        stderr = zeros(nGroups, numel(bins));
    end
    for g = 1:nGroups
        this_hist = baSm(myData.trGroups==trGroupLabels(g),:);
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
else
    psthSm = mean(baSm);
    if myData.params.showErrorShading
        stderr = std(baSm)./sqrt(size(baSm,1));
    end
end

% compute raster
if myData.params.showAllTraces
    [~, inds] = sort(myData.trGroups);
    [tr,b] = find(ba(inds,:));
else
    [tr,b] = find(ba);
end
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything

% scale the raster ticks
rasterY(2:3:end) = rasterY(2:3:end)+myData.params.rasterScale;

% compute the tuning curve
tuningCurve = zeros(nGroups,2);
for g = 1:nGroups
    theseCounts = spikeCounts(myData.trGroups==trGroupLabels(g));
    tuningCurve(g,1) = mean(theseCounts);
    tuningCurve(g,2) = std(theseCounts)./sqrt(length(theseCounts));
end

%% ORGANIZE BEHAVIORAL TRIAL DATA & COMPUTE PSTH

if strcmp(myData.expt_type,'free')
    myData.trialTime = myData.sync_data.trial_on(myData.params.trialIndex);
    myData.trialLength = myData.sync_data.trial_off(myData.params.trialIndex)-myData.trialTime;
    
    % get firing rate trace for entire trial
    iTrial = myData.params.trialIndex;
    trial_bin_size = .1;
    trial_spikes = st(st>myData.trialTime & st<(myData.trialTime+myData.trialLength));
    [trial_ba, trial_bin_centers] = timestampsToBinned(trial_spikes,myData.trialTime,trial_bin_size, [0 myData.trialLength]);
    trial_psth = trial_ba./trial_bin_size; % normalize to Hz
    % PSTH smoothing filter
    trial_gw = gausswin(8);
    trial_smWin = trial_gw./sum(trial_gw);
    trial_baSm = conv2(trial_smWin,1,trial_psth', 'same')'./trial_bin_size; % smooth ba
    
    % get behavioral data from etho
    etho = myData.etho;
    trial_time = etho.trial_time{iTrial};
    velocity = etho.velocity{iTrial};
    idx = isnan(velocity);
    velocity(idx) = 0; % i'm setting these to 0 right now as a quick fix.
    
    x_center = etho.x_center{iTrial};
    y_center = etho.y_center{iTrial};
    pd_stamps = trial_time(find(etho.etho_photodiode{iTrial}));
    nest_stamps = trial_time(find(etho.in_nest{iTrial}));
    grooming_stamps = trial_time(find(etho.grooming{iTrial}));
    rearing_stamps = trial_time(find(etho.rearing{iTrial}));
    eating_stamps = trial_time(find(etho.eating{iTrial})); %#ok<*FNDSB>
    bout_stamps = trial_time(find(etho.bout{iTrial}));
    exploring_stamps = trial_time(find(etho.exploring{iTrial}));
end

%% EXPT INFORMATION
expt_string = sprintf('%d total neurons',length(myData.clusterDepths));
if median(myData.sync_data.trialPulses) == 3
    trial_string = sprintf('%d trials',length(myData.sync_data.stimIDs)/3);
else
    trial_string = sprintf('%d trials',length(myData.sync_data.stimIDs));
end
exptInfo = sprintf('%s\n\n%s\n\n%s',myData.meta.name,expt_string,trial_string);

% plot in figure
expt_panel = myData.panels(4);
expt_text = uicontrol('Style','text','String',exptInfo,'parent',expt_panel,'Position',[10 40 120 120],'BackgroundColor','w','FontSize',10,'HorizontalAlignment','Left');


%% CLUSTER INFORMATION
iClu = myData.params.clusterIndex;
depth = myData.clusterDepths(iClu);
group = myData.cgs(iClu);
this_cid = myData.cids(iClu);

if group == 2
    group = 'Good';
elseif group == 1
    group = 'MUA';
end

% these borders are for NP8!
v1Borders = myData.borders(1,:);
hcBorders = myData.borders(2,:);

if size(myData.borders,1) >= 3
    mbBorders = myData.borders(3,:);
else
    mbBorders = [];
end


% would be nice to have brain region here too, but it needs to be in the
% meta data
depth_string = sprintf('Distance from tip: %2.2f mm',(depth/1000));
id_string = sprintf('Cluster ID: %d',this_cid);
txtInfo = sprintf('%s\n\n%s\n%s',depth_string,group,id_string);

% plot in figure
id_panel = myData.panels(1);
id_text = uicontrol('Style','text','String',txtInfo,'parent',id_panel,'Position',[10 40 100 120],'BackgroundColor','w','FontSize',10, 'HorizontalAlignment','left');

%% GET WAVEFORM
% for now, i'm going to use the template from kilosort.
% importing raw WFs is much more memory consuming...
% medWFs = extractMedianWFs(clu, st, Fs, data_path, 'int16', [nCh nSamps], chanMap, [1,1])
% spike_templates gives the template for each spike
% clu gives the clu for each spike

[cluIDs] = find(myData.clu == this_cid);
templateID = myData.spikeTemplates(cluIDs(1))+1; % it doesn't actually matter which cluID(i) we use, they should be the same
wf_panel = myData.panels(2);
axes(myData.plotAxes(4));
hold off;
plot(myData.waveforms(templateID,:));
set(gca,'YTickLabel',[],'XTickLabel',[],'XLim',[0 82],'Visible','off');
axis tight;

%% PLOT

% THIS IS WHERE YOU DETERMINE THE COLOR OF EACH PLOT
switch myData.stim_type
    case 'AVLR'
        if strcmp(myData.expt_type, 'free')
            if length(unique(myData.trGroups)) == 3
                colors(1,:) = [.8 .3 .1]; % high AV
                colors(2,:) = [0 .3 .4]; % high visual
                colors(3,:) = [0 0 0]; % sound only
            else
                colors(1,:) = [.8 .3 .1]; % high AV
                colors(2,:) = [0 .3 .4]; % high visual
                colors(3,:) = [.9 .6 .3]; % low AV
                colors(4,:) = [0 .5 .8]; % low visual
                colors(5,:) = [0 0 0]; % sound only
                colors(6,:) = [.5 .5 .5]; % blank
            end
        else
            colors = parula(nGroups);
        end
    case 'ori'
        colors = parula(nGroups);
end

%% PSTH
axes(myData.plotAxes(1));
hold off;
if myData.params.showAllTraces
    for g = 1:nGroups
        % plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
        psth_line = shadedErrorBar(bins,psthSm(g,:),stderr(g,:));
        psth_line.mainLine.Color = colors(g,:);
        psth_line.patch.FaceColor = colors(g,:);
        psth_line.patch.FaceAlpha = .3;
        psth_line.patch.LineStyle = 'none';
        psth_line.mainLine.LineWidth = 2;
        hold on;
    end
else
    %plot(bins, psthSm);
    psth_line = shadedErrorBar(bins,psthSm,stderr);
    psth_line.patch.FaceAlpha = .3;
    psth_line.patch.LineStyle = 'none';
    psth_line.mainLine.LineWidth = 1;
    % ylim([0 38])
    % xticks([-.5:.1:2])
    
end
xlim(myData.params.window);
title(['cluster ' num2str(myData.clusterIDs(myData.params.clusterIndex))]);
xlabel('time (sec)');
ylabel('firing rate (Hz)');


% some day we'll figure out this legend situation
% legend_text = sprintf('%s\n%s','high AV','high vis');
% legend = uicontrol('Style','text','Parent',myData.panels(5),'String',legend_text,'Position',[50 50 100 120],'BackgroundColor','transparent','FontSize',10, 'HorizontalAlignment','left');
% legend({'high AV','high vis','low AV','low vis','sound','blank'})
% legend('Location','best'); legend(gca,'boxoff')
set(gca,'FontSize',8)
box off;

%% RASTER
axes(myData.plotAxes(2));
hold off;
plot(rasterX,rasterY, 'k');
xlim(myData.params.window);
ylim([0 length(myData.eventTimes)+1]);
ylabel('event number'); xlabel('time (sec)');
set(gca,'FontSize',8)
set(get(gca, 'Title'), 'FontSize', 8);
box off;
makepretty;

%% TUNING CURVE
axes(myData.plotAxes(3));
hold off;
switch myData.stim_type
    case 'AVLR'
        AVLRbar = barwitherr(tuningCurve(:,2),tuningCurve(:,1));
        if length(unique(myData.trGroups)) == 3
            set(gca,'XTickLabel',{'high AV','high vis','sound'},'FontSize',8)
        else
            set(gca,'XTickLabel',{'high AV','high vis','low AV','low vis','sound','blank'},'FontSize',8)
        end
        set(AVLRbar(1),'FaceColor',[.5 .5 .5]) % would be nice if these colors matched AVLR colors
        xlabel('stimulus type')
    case 'ori'
        errorbar(trGroupLabels, tuningCurve(:,1), tuningCurve(:,2), 'o-');
        xlabel('orientation');
end
ylabel('avg firing rate (Hz)');
set(gca,'FontSize',8)
box off;

%% TRIAL FIRING RATE

if strcmp(myData.expt_type,'free')
    axes(myData.plotAxes(7));
    hold off;
    plot(trial_baSm,'LineWidth',2);
    set(gca,'XTick',[],'FontSize',8)
    title(['trial ' num2str(myData.params.trialIndex)]);
    xlabel('time (same as above)');
    ylabel('firing rate (Hz) (.1 ms bins)');
    box off; axis tight;
    
    %% VELOCITY
    axes(myData.plotAxes(6));
    hold off;
    plot(velocity,'LineWidth',2,'Color',[.8 .2 0]);
    
    % figure out trial type
    switch myData.sync_data.stimIDs(myData.params.trialIndex)
        case 1; trial_type = 'high AV';
        case 2; trial_type = 'high visual';
        case 3; trial_type = 'low AV';
        case 4; trial_type = 'low visual';
        case 5; trial_type = 'sound only';
        case 6; trial_type = 'blank';
    end
    title(['trial type:' trial_type]);
    set(gca,'XTick',[],'FontSize',8);box off; axis tight; ylabel('Velocity (cm/s)')
    
    %% TRIAL BEHAVIOR -- plot available behavioral variables
    axes(myData.plotAxes(5));
    cla reset
    hold off;
    if ~isempty(eating_stamps)
        eating_patch = patch([eating_stamps(1),eating_stamps(1),eating_stamps(end),eating_stamps(end)],[0,max(velocity),max(velocity),0],[0.2656    0.7148    0.1406]);
        set(eating_patch,'LineStyle','none');
    else
        eating_patch = patch([0,0,1,1],[0,1,1,0],[0.2656    0.7148    0.1406]);
        set(eating_patch,'FaceAlpha',0,'LineStyle','none')
    end
    hold on;
    if ~isempty(bout_stamps)
        bout_patch = patch([bout_stamps(1),bout_stamps(1),bout_stamps(end),bout_stamps(end)],[0,max(velocity),max(velocity),0],[0.2578    0.5234    0.9531]);
        set(bout_patch,'LineStyle','none')
    else
        bout_patch = patch([0,0,1,1],[0,1,1,0],[0.2578    0.5234    0.9531]);
        set(bout_patch,'FaceAlpha',0,'LineStyle','none')
    end
    hold on;
    if ~isempty(rearing_stamps)
        rearing_patch = patch([rearing_stamps(1),rearing_stamps(1),rearing_stamps(end),rearing_stamps(end)],[0,max(velocity),max(velocity),0],[0.9258    0.7969    0.1641]);
        set(rearing_patch,'LineStyle','none')
    else
        rearing_patch = patch([0,0,1,1],[0,1,1,0],[0.9258    0.7969    0.1641]);
        set(rearing_patch,'FaceAlpha',0,'LineStyle','none')
    end
    if ~isempty(nest_stamps)
        nest_patch = patch([nest_stamps(1),nest_stamps(1),nest_stamps(end),nest_stamps(end)],[0,max(velocity),max(velocity),0],[.8 .8 .8]);
        set(nest_patch,'LineStyle','none')
    else
        nest_patch = patch([0,0,1,1],[0,1,1,0],[.8 .8 .8]);
        set(nest_patch,'FaceAlpha',0,'LineStyle','none')
    end
    if ~isempty(exploring_stamps)
        exploring_patch = patch([exploring_stamps(1),exploring_stamps(1),exploring_stamps(end),exploring_stamps(end)],[0,max(velocity),max(velocity),0],[0.8750    0.5391    0.1016]);
        set(exploring_patch,'LineStyle','none')
    else
        exploring_patch = patch([0,0,1,1],[0,1,1,0],[0.8750    0.5391    0.1016]);
        set(exploring_patch,'FaceAlpha',0,'LineStyle','none')
    end
    if ~isempty(pd_stamps)
        pd_patch = patch([pd_stamps(1),pd_stamps(1),pd_stamps(end),pd_stamps(end)],[0,max(velocity),max(velocity),0],[0 0 0]);
        set(nest_patch,'LineStyle','none')
    else
        pd_patch = patch([0,0,1,1],[0,1,1,0],[0 0 0]);
        set(pd_patch,'FaceAlpha',0,'LineStyle','none')
    end
    
    hold off; box off;
    xlabel('time (secs)'); set(gca,'YTick',[])
    legend({'eating','bout','rearing','nest','exploring','photodiode'})
    %legend(gca,'boxoff')
    legend(gca,'Location','best')
    axis tight
end

end

function psthViewerCallback(f, keydata)

updateOtherFigs = false;

myData = get(f, 'UserData');

switch keydata.Key
    case 'rightarrow' % increment cluster index
        
        myData.params.clusterIndex = myData.params.clusterIndex+1;
        if myData.params.clusterIndex>length(myData.clusterIDs)
            myData.params.clusterIndex=1;
        end
        updateOtherFigs = true;
        
    case 'leftarrow' % decrement cluster index
        
        myData.params.clusterIndex = myData.params.clusterIndex-1;
        if myData.params.clusterIndex<1
            myData.params.clusterIndex=length(myData.clusterIDs);
        end
        updateOtherFigs = true;
        
    case 'uparrow' % increase smoothing
        myData.params.smoothSize = myData.params.smoothSize*1.2;
        
    case 'downarrow' % decrease smoothing
        myData.params.smoothSize = myData.params.smoothSize/1.2;
        
    case 'e' % whether to show standard error as shading
        myData.params.showErrorShading = ~myData.params.showErrorShading;
        
    case 'p' % whether to plot the psth trace for each condition or just the overall one
        myData.params.showAllTraces = ~myData.params.showAllTraces;
        updateOtherFigs = true;
        
    case 'r'
        ax = myData.plotAxes(2);
        title('click start and stop of range');
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        myData.params.startRange = q(1,1);
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        myData.params.stopRange = q(1,1);
        if myData.params.stopRange<myData.params.startRange
            tmp = myData.params.startRange;
            myData.params.startRange = myData.params.stopRange;
            myData.params.stopRange = tmp;
        end
        
    case 'c'
        newC = inputdlg('cluster ID?');
        ind = find(myData.clusterIDs==str2num(newC{1}),1);
        if ~isempty(ind)
            myData.params.clusterIndex = ind;
        end
        
        updateOtherFigs = true;
        
    case 't'
        newT = inputdlg('trial #?');
        myData.params.trialIndex = str2num(newT{1});
        
        updateOtherFigs = true;
        
end

set(f, 'UserData', myData);

% plot with new settings
psthViewerPlot(f)

if updateOtherFigs && f==keydata.Source.Number && isfield(myData, 'otherFigs')
    % checking that the current figure matches the source number prevents
    % doing this when called *not* as the original fig
    setOtherFigsClusterIndex(myData, myData.params.clusterIndex)
    plotOtherFigs(f)
end

end


function setOtherFigsClusterIndex(myData, cInd)

for thatFig = myData.otherFigs
    thatData = get(thatFig, 'UserData');
    thatData.params.clusterIndex = cInd;
    set(thatFig, 'UserData', thatData);
    
end

end

function plotOtherFigs(f)
myData = get(f, 'UserData');
for thatFig = myData.otherFigs
    %     thatFigFcn = get(thatFig, 'KeyPressFcn');
    figs = get(0, 'Children');
    figNumsCell = get(figs, 'Number');
    figNums = [figNumsCell{:}];
    thatFigObj = figs(figNums==thatFig);
    psthViewerPlot(thatFigObj);
end
figure(f) % return focus here
end