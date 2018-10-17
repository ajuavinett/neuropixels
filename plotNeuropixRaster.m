function plotNeuropixRaster(meta,sp,win,sync_data)

% to plot raster of neuropix data with designated size
% window and raster size
% e.g. win = [1306 1308];

%% SCRIPT OPTIONS
rasterScale = 20;
border_bit = 1;
save_bit = 0;
sync_bit = 1;

%% GET INFO OUT OF SP STRUCTURE

% f = figure; set(f, 'Color', 'w');
st = sp.st;
inclSpikes = st>win(1) & st<=win(2);
st = st(inclSpikes);
clu = sp.clu;
clu = clu(inclSpikes);
cids = sp.cids;
cgs = sp.cgs;
sd = sp.spikeDepths;
sd = sd(inclSpikes);

co = get(gca, 'ColorOrder'); nc = size(co,1);

%% PLOT

for c = 1:length(sp.cids)
    if sp.cgs(c)==2
        thisColor = co(mod(c,nc)+1,:);
    else
        thisColor = [0.7 0.7 0.7]; % grey for unsorted and mua
    end
    
    [xx, yy] = rasterize(st(clu==cids(c)));
    theseSD = sd(clu==cids(c));
    yy(1:3:end) = yy(1:3:end)+theseSD';
    yy(2:3:end) = yy(2:3:end)*rasterScale+theseSD';
    
    plot(xx,yy, 'Color', thisColor, 'LineWidth', 1.5);
    hold on;
end

box off;

if border_bit
    for iBorder = 1:(size(meta.borders,2)-1)
        plot([win(1) win(end)], meta.borders(iBorder,1)*[1 1], 'k--', 'LineWidth', 1.0);
    end
end

%% added this so you can add timestamps for when stimuli play
if sync_bit
    syncIDs = find(sync_data.photodiode>win(1) & sync_data.photodiode<win(2));
    for iSync = 1:length(syncIDs)
        plot(sync_data.photodiode(syncIDs(iSync))*[1 1],[0 4000],'k','LineWidth',1)
    end
end

set(gca, 'XTick', []);
ylim([0 max(sp.ycoords)+20])
ylabel('depth on probe (um)')
xlabel('time (sec)');

%% SAVE FIGURE
% if save_bit
%     savedir = fullfile(meta.datadir,'Figures');
%
%     if ~exist(savedir,'dir')
%         mkdir(savedir)
%     end
%
%     % plot2svg(fullfile(datadir,'neuronsbydepth'),gcf)
%     saveas(gcf,strcat(savedir,'\','raster'),'fig')
%     saveas(gcf,strcat(savedir,'\','raster'),'svg')
% end
