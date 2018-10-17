function [areaCounts,clusterDepths_good] = plotClustersbyDepth(meta,sp)

% basic plot of clusters over depth
% on each probe, higher depth numbers are superficial, i.e. nearer to the
% top of the brain; lower numbers are deeper, nearer the tip

% f = figure;
% set(gcf,'Position',[680   558   136   420],'color','w')

clusterDepths = sp.clusterDepths;
ca = sp.clusterAmps;
cgs = sp.cgs;

xx = rand(size(cgs));

scatter(xx(cgs==1),clusterDepths(cgs==1),ca(cgs==1)/5); % plot MUA
hold on;
scatter(xx(cgs==2),clusterDepths(cgs==2),ca(cgs==2)/5) % plot good units
title(meta.name,'interpreter','none');
set(gca,'XTickLabel',[])
ylabel('depth on probe (um)')
% legend({'MUA', 'Good'});
ylim([0 max(sp.ycoords)+20])
box off


nBorders = size(meta.borders,1);
% meta.borders = meta.borders{1};

for iBorder = 1:nBorders
    plot([0 1], meta.borders(iBorder,1)*[1 1], 'k--', 'LineWidth', 1.0);
end

% get area counts based on given borders
if nBorders >= 3
    areaCounts(1) = sum(clusterDepths>=meta.borders(1,1) & cgs==2);
    areaCounts(2) = sum(clusterDepths>=meta.borders(2,1) & clusterDepths<meta.borders(1,1) & cgs==2);
    areaCounts(3)= sum(clusterDepths<meta.borders(2,1) & cgs==2);
    fprintf(1, 'recorded %d v1 neurons, %d hippocampus, %d midbrain\n', areaCounts(1), areaCounts(2), areaCounts(3));
elseif nBorders == 2
    areaCounts(1) = sum(clusterDepths>=meta.borders(1,1) & cgs==2);
    areaCounts(2)= sum(clusterDepths<meta.borders(1,1) & cgs==2);
    fprintf(1, 'recorded %d SC neurons, %d midbrain\n', areaCounts(1), areaCounts(2));
end

clusterDepths_good = clusterDepths(cgs==2); % only save good clusters


%% SAVE FIGURE

% savedir = fullfile(meta.datadir,'Figures');
%
% if ~exist(savedir,'dir')
%     mkdir(savedir)
% end

% plot2svg(fullfile(datadir,'neuronsbydepth'),gcf)
% saveas(gcf,strcat(savedir,'\','neuronsbydepth'),'fig')
% saveas(gcf,strcat(savedir,'\','neuronsbydepth'),'svg')