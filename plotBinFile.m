function plotBinFile(filename,opt)
% filename = 'C:\Neuropixels\NP6\run3_g0_t0.imec.ap.bin';
% to quickly plot a chunk of a bin file and look at sample traces

%% FILE OPTIONS
plot_all_bit = 0; % imagesc plot
separate_bit = 1; % chooses random traces to plot separately
overlap_bit = 0; % plots those same random traces but overlapped
start_time = 25000000;
end_time = 25001000;
chunkSize = 100000000;

if opt == 4
    nChan = 277;
else
    nChan = 385; %384 + sync?!
end

chosenCh = [103,139,140,141]; % choose 4-6 random channels to plot

%% READ FILE
if nargin < 1
    [FileName,PathName,FilterIndex] = uigetfile('*.bin','Pick a bin file');
    filename = fullfile(PathName, FileName);
    if ~FilterIndex
        return;
    end
end

fid = fopen(filename, 'r');

if fid == 3
    disp('successfully opened bin file')
end

dat = fread(fid, [nChan chunkSize], '*int16');

fclose(fid);

% chanMap = readNPY('channel_map.npy');
% dat = dat(chanMap+1,:);

%% PLOT

% plot all
if plot_all_bit
    figure; imagesc(dat(:,1:30000))
end

% plot some random traces
if separate_bit
    figure;
    set(gcf,'Position',[6 85 1908 893]);

    for iCh = 1:length(chosenCh)
        subplot(length(chosenCh),1,iCh);
        % scrollplot(plot(dat(chosenCh(iCh),[start_time:end_time])),'WindowSizeX',1000)
        plot(dat(chosenCh(iCh),[start_time:end_time]))
        xlabel('time')
        ylabel(num2str(chosenCh(iCh)))
    end
end

% to plot channels overlapped (this is good to check that theyre not offset
% in time
if overlap_bit
    figure;
    for iCh = 1:length(chosenCh)
        scrollplot(plot(dat(chosenCh(iCh),[start_time:end_time])),'WindowSizeX',1000)
        hold on
        xlabel('time')
        ylabel(num2str(chosenCh(iCh)))
    end
end