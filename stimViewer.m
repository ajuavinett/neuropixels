function stimViewer(sp,sync_data,meta)
%
% Controls:
% - c: dialog box to pick a new cluster ID number
% - t: toggle showing psth traces for each grouping variable or just the
% overall. If showing just overall, raster is sorted chronologically. If
% showing by grouping variable, raster is sorted by that variable.
% - r: select a new range within which to count spikes for the tuning curve

%% SET UP FIGURE
f = figure;
set(gcf,'Position',[1     39    1619    858],'Color','w')
expt_panel = uipanel('Title','Expt Info','FontSize',12,'BackgroundColor','white',...
    'Position',[.05 .65 .3 .3]);
id_panel = uipanel('Title','Cluster Info','FontSize',12,'BackgroundColor','white',...
    'Position',[.05 .35 .3 .3]);
wf_panel= uipanel('Parent',id_panel,'Title','Template','FontSize',12,...
    'BackgroundColor','w','Position',[.5 .1 .4 .9]);
wf_ax = axes('parent',wf_panel);
psth_panel = uipanel('Title','PSTH for all groups','FontSize',12,'BackgroundColor','white',...
    'Position',[.35 .65 .3 .3]);
psth_ax = axes('parent',psth_panel);
raster_panel = uipanel('Title','Raster Plot','FontSize',12,'BackgroundColor','white',...
    'Position',[.35 .35 .3 .3]);
raster_ax = axes('parent',raster_panel);
tuning_panel = uipanel('Title','Tuning','FontSize',12,'BackgroundColor','w',...
    'Position',[.35 .05 .3 .3]);
tuning_ax = axes('parent',tuning_panel);
velocity_panel = uipanel('Title','Velocity','FontSize',12,'BackgroundColor','white',...
    'Position',[.05 .05 .3 .3]);
velocity_ax = axes('parent',velocity_panel);
type_panel = uipanel('Title','Grouped by Looming vs Receding','FontSize',12,'BackgroundColor','white',...
    'Position',[.65 .05 .3 .3]);
type_ax = axes('parent',type_panel);
modality_panel = uipanel('Title','Grouped by Modality','FontSize',12,'BackgroundColor','white',...
    'Position',[.65 .35 .3 .3]);
modality_ax = axes('parent',modality_panel);
contrast_panel = uipanel('Title','Grouped by Contrast','FontSize',12,'BackgroundColor','white',...
    'Position',[.65 .65 .3 .3]);
contrast_ax = axes('parent',contrast_panel);

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
myData.sp = sp;

eventTimes = sync_data.photodiode;
if ~issorted(eventTimes)
    [eventTimes, ii] = sort(eventTimes);
    trGroups = trGroups(ii);
end

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
myData.plotAxes = [psth_ax,raster_ax,tuning_ax,wf_ax,velocity_ax,contrast_ax,modality_ax,type_ax];
myData.panels = [id_panel,wf_panel,expt_panel,psth_panel,contrast_panel,modality_panel,type_panel];
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

%% GET PSTH FOR EACH TRIAL TYPE
if myData.params.showAllTraces
    [psthSm, stderr] =  generatePSTHforGroup(baSm,bins,myData,myData.trGroups,trGroupLabels);
    if length(unique(myData.trGroups)) == 9 %NP16
        [psthSm_type,stderr_type] = generatePSTHforGroup(baSm,bins,myData,myData.sync_data.sound_dir,[2 3]);    
    else
        [psthSm_type,stderr_type] = generatePSTHforGroup(baSm,bins,myData,myData.sync_data.type,[1 2]);
    end
    [psthSm_contrast,stderr_contrast] = generatePSTHforGroup(baSm,bins,myData,myData.sync_data.contrast,unique(myData.sync_data.contrast));
    [psthSm_modality,stderr_modality] = generatePSTHforGroup(baSm,bins,myData,myData.sync_data.sound_bit,[0 1]); % JUST SOUNDBIT, CURRENTLY
else
    psthSm = mean(baSm);
    if myData.params.showErrorShading
        stderr = std(baSm)./sqrt(size(baSm,1));
    end
end

% compute raster
if myData.params.showAllTraces
    [~, inds] = sort(myData.trGroups); % raster is sorted by groups
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

%% EXPT INFORMATION PANEL
expt_string = sprintf('%d total neurons',length(myData.clusterDepths));
if median(myData.sync_data.trialPulses) == 3
    trial_string = sprintf('%d trials',length(myData.sync_data.stimIDs)/3);
else
    trial_string = sprintf('%d trials',length(myData.sync_data.stimIDs));
end
exptInfo = sprintf('%s\n\n%s\n\n%s',myData.meta.name,expt_string,trial_string);

expt_panel = myData.panels(3);
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

%% PSTH FOR ALL CURVES
axes(myData.plotAxes(1));
if myData.nGroups == 6
    groupNames = ['HI AVL','HI VL ','SOUND ','BLANK ','HI AVR','HI VR '];
elseif myData.nGroups == 8
    groupNames = ['HI AVL','HI VL ','LOWAVL','LOW VL','HI AVR','HI VR ','LOWAVR','LOW VR'];
elseif myData.nGroups == 14 %NP9
    groupNames = ['HI AVL','HI VL ','LOWAVL','LOW VL','SOUND ','BLANK ','HI AVR','HI VR ','LOWAVR','LOW VR','MID VL','MID VR','MIDAVL','MIDAVR'];
elseif myData.nGroups == 9 %NP16
    groupNames = ['HIVLAL','HIVLAR','HI VL ','AL    ','AR    ','MIVLAL','MIVLAR','MID VL','BLANK '];
else
    errordlg('please define groups')
end

hold off;
if myData.params.showAllTraces
    for g = 1:nGroups
        colors = parula(nGroups);
        % plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
        psth_line = shadedErrorBar(bins,psthSm(g,:),stderr(g,:));
        psth_line.mainLine.Color = colors(g,:);
        psth_line.patch.FaceColor = colors(g,:);
        psth_line.patch.FaceAlpha = .3;
        psth_line.patch.LineStyle = 'none';
        psth_line.mainLine.LineWidth = 2;
        psth_line.mainLine.DisplayName = groupNames((g-1)*6+1:6*g);
        legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','EastOutside');
        legend('boxoff')
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
set(gca,'FontSize',8)
box off;

%% PSTH FOR TYPE
axes(myData.plotAxes(8));
if length(unique(myData.trGroups)) == 9 %NP16
    groupNames = ['LOOMING  AUD','RECEDING AUD'];
else
    groupNames = ['LOOMING  VIS','RECEDING VIS'];
end
hold off;
for g = 1:2
    type_colors = [0,0,0;.8,.8,.8];
    psth_line = shadedErrorBar(bins,psthSm_type(g,:),stderr_type(g,:));
    psth_line.mainLine.Color = type_colors(g,:);
    psth_line.patch.FaceColor = type_colors(g,:);
    psth_line.patch.FaceAlpha = .3;
    psth_line.patch.LineStyle = 'none';
    psth_line.mainLine.LineWidth = 2;
    psth_line.mainLine.DisplayName = groupNames((g-1)*12+1:12*g);
    legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','EastOutside');
    legend('boxoff')
    hold on;
end

xlim(myData.params.window);
xlabel('time (sec)');
ylabel('firing rate (Hz)');
set(gca,'FontSize',8)
box off;

%% PSTH FOR SOUND BIT
axes(myData.plotAxes(7));
groupNames = ['VISUAL ONLY','AUDIOVISUAL'];
hold off;
for g = 1:2
    type_colors = [.8,.8,.8;0,0,0];
    psth_line = shadedErrorBar(bins,psthSm_modality(g,:),stderr_modality(g,:));
    psth_line.mainLine.Color = type_colors(g,:);
    psth_line.patch.FaceColor = type_colors(g,:);
    psth_line.patch.FaceAlpha = .3;
    psth_line.patch.LineStyle = 'none';
    psth_line.mainLine.DisplayName = groupNames((g-1)*11+1:11*g);
    legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','EastOutside');
    legend('boxoff')
    psth_line.mainLine.LineWidth = 2;
    hold on;
end

xlim(myData.params.window);
xlabel('time (sec)');
ylabel('firing rate (Hz)');
set(gca,'FontSize',8)
box off;

%% PSTH FOR CONTRAST
axes(myData.plotAxes(6));
if length(unique(myData.sync_data.contrast)) > 1 % if we tested more than one contrast
    
    if length(unique(myData.sync_data.contrast)) == 2
        groupNames = ['LOW  ','HIGH '];
    elseif length(unique(myData.sync_data.contrast)) == 3
        groupNames = ['NONE ','LOW  ','HIGH '];
    elseif length(unique(myData.sync_data.contrast)) == 4 %NP9
        groupNames = ['NONE ','LOW  ','MID  ','HIGH '];
    end
    
    contrast_colors = flip(gray(size(stderr_contrast,1)));
    hold off;
    for g = 1:length(unique(myData.sync_data.contrast))
        psth_line = shadedErrorBar(bins,psthSm_contrast(g,:),stderr_contrast(g,:));
        psth_line.mainLine.Color = contrast_colors(g,:);
        psth_line.patch.FaceColor = contrast_colors(g,:);
        psth_line.patch.FaceAlpha = .3;
        psth_line.patch.LineStyle = 'none';
        psth_line.mainLine.LineWidth = 2;
        psth_line.mainLine.DisplayName = groupNames((g-1)*5+1:5*g);
        legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','EastOutside');
        legend('boxoff')
        hold on;
    end
    
    xlim(myData.params.window);
    xlabel('time (sec)');
    ylabel('firing rate (Hz)');
    set(gca,'FontSize',8)
    box off;
end

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
        elseif length(unique(myData.trGroups)) == 6 && max(unique(myData.trGroups)) == 8
            set(gca,'XTickLabel',{'high AV looming','high visual looming','sound','blank','high AV receding','high visual receding'},'FontSize',8)
        elseif length(unique(myData.trGroups)) == 8
            set(gca,'XTickLabel',{'high AV looming','high visual looming','low AV looming','low visual looming ','high AV receding','high visual receding','low AV receding','low visual receding'})
        elseif length(unique(myData.trGroups)) == 14
            set(gca,'XTickLabel',{'HI AVL','HI VL ','LOWAVL','LOW VL','SOUND ','BLANK ','HI AVR','HI VR ','LOWAVR','LOW VR','MID VL','MID VR','MIDAVL','MIDAVR'});
        elseif length(unique(myData.trGroups)) == 9
            set(gca,'XTickLabel',{'HIVLAL','HIVLAR','HI VL ','AL    ','AR    ','MIVLAL','MIVLAR','MID VL','BLANK '});
        else
            errordlg('need to fix labels for tuning curve plot')
        end
        set(AVLRbar(1),'FaceColor',[.5 .5 .5]) % would be nice if these colors matched AVLR colors
        xlabel('stimulus type')
    case 'ori'
        errorbar(trGroupLabels, tuningCurve(:,1), tuningCurve(:,2), 'o-');
        xlabel('orientation');
end
ylabel('avg firing rate (Hz)');
xtickangle(30);
set(gca,'FontSize',8)
box off;

%% WHEEL RUNNING SPEED VS VELOCITY

[binnedRates,velocity,~] = mvmtNeurons(myData.sp,myData.sync_data,0);
velocity_scale_factor = 100;
velocity_window = [1:200];
axes(myData.plotAxes(5));
hold off
iNeuron = myData.params.clusterIndex;
nBinnedRatesbins = length(binnedRates);
binnedRates = binnedRates(:,[1:nBinnedRatesbins-1]);
plot(velocity_window,velocity(velocity_window)*velocity_scale_factor,'LineWidth',2,'Color',[.8 .8 .8])
hold on
plot(velocity_window,binnedRates(iNeuron,velocity_window),'LineWidth',2)
hold off

xlabel('time (samples)');
ylabel('firing rate (Hz) || running speed');
box off; axis tight;

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