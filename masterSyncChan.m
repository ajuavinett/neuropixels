function sync_data = masterSyncChan(opt,folder)


%% CHOOSE SYNC FILE
% note, it is much faster if the LFP file is local
if nargin < 2
    folder = uigetdir('C:\','Choose directory of LFP file');
end

%% get info about experiment based on probe option

if opt == 4
    numChans = 277;
else
    numChans = 385;
end

Fs = 2500; % sampling rate of LFP file, which we'll use to extract sync
syncChanIndex = numChans; % it should always be the last channel

%% extract sync channel
sync_data = {};
syncDat = extractSyncChannel(folder, numChans, syncChanIndex);

% get times in seconds. eventTimes{1} = [on,off], {2} = on, {3} = off
[eventTimes, bitRep] = spikeGLXdigitalParse(syncDat, Fs);


%% hard coded for which pins these are on basestation sync board
trial_pin = 1;
sync_pin = 2;
photodiode_pin = 3;
center_pin = 7;
nest_pin = 5;
mvmt_pin1 = 9;
mvmt_pin2 = 11;

%% save on/off times and traces

% get trial start & end times
trial_on = eventTimes{trial_pin}(2);
trial_off = eventTimes{trial_pin}(3);

% recreate sync trace
sync_on = eventTimes{sync_pin}(1); % get sync on
sync_trace = bitRep(:,sync_pin); % save entire trace
% figure; set(gcf,'Position',[8 558 1908 420]); % to plot sync trace
% scrollplot(plot(sync_trace),'WindowSize','1000');

%% get photodiode, center, nest times
photodiode = eventTimes{photodiode_pin}(3); % 3 is the "on" signal (the PD trace goes DOWN)
%photodiode_on = eventTimes{photodiode_pin}(2);
%photodiode_off = eventTimes{photodiode_pin}(3);
photodiode_trace = bitRep(:,photodiode_pin); % save entire trace
center = eventTimes{center_pin}(2);
nest = eventTimes{nest_pin}(2);

%% get wheel mvmt times

mvmt1_trace = bitRep(:,mvmt_pin1); % save entire trace 1
mvmt2_trace = bitRep(:,mvmt_pin2); % save entire trace 2

if min(mvmt1_trace) == max(mvmt1_trace) % for later experiments, the pin is different
    mvmt_pin1 = 6;
    mvmt_pin2 = 7;
    mvmt1_trace = bitRep(:,mvmt_pin1); % save entire trace 1
    mvmt2_trace = bitRep(:,mvmt_pin2); % save entire trace 2    
end

mvmt1_on = eventTimes{mvmt_pin1}(2);
mvmt1_on = mvmt1_on{1}*Fs*1000 ; % get it back to samples
mvmt1_off = eventTimes{mvmt_pin1}(3);
mvmt1_off = mvmt1_off{1}*2500*1000; % get it back to samples
mvmt2_on = eventTimes{mvmt_pin2}(2);
mvmt2_on = mvmt2_on{1}*Fs*1000 ; % get it back to samples
mvmt2_off = eventTimes{mvmt_pin2}(3);
mvmt2_off = mvmt2_off{1}*2500*1000; % get it back to samples



%% save as sync_data struct
sync_data.trial_on = trial_on{1};
sync_data.trial_off = trial_off{1};
sync_data.sync_on = sync_on;
sync_data.sync_trace = sync_trace;
sync_data.photodiode = photodiode{1};
sync_data.pd_trace = photodiode_trace;
sync_data.center = center{1};
sync_data.nest = nest{1};
sync_data.mvmt1_on = mvmt1_on;
sync_data.mvmt1_trace = mvmt1_trace;
sync_data.mvmt2_on = mvmt2_on;
sync_data.mvmt2_trace = mvmt2_trace;

%% CLEAN UP SYNC FILE
% for some reason, sometimes there are syncs in the very beginning
% (and this is impossible)

if exist('sync_data.trial_on','var') && round(sync_data.trial_on(1)) == 0
    sync_data.trial_on = sync_data.trial_on(2:end);
    sync_data.trial_off = sync_data.trial_off(2:end);
end

if round(sync_data.photodiode(1)) == 0
    sync_data.photodiode = sync_data.photodiode(2:end);
end

%% SAVE SYNC FILE

save('sync_data.mat','sync_data')



