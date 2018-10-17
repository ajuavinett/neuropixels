function edited_sync_data = makeHeadfixedSync(expt_name,sync_data,expt_type)

% THIS SCRIPT IN PROGRESS
% expt_type = 'ori' or 'retinotopy' or 'CSD'
% sync_raw is the output of masterSyncChan (which you should have run)

photodiode = sync_data.photodiode;
[expt_parts,~] = strsplit(expt_name,'_');
anim = expt_parts{1};
expt = expt_parts{2};

% import analyzer
[trialArray,headers,~] = makeTrialArray(anim,expt);

switch expt_type
    
    case 'AVLR'
        for iHeader = 1:length(headers)
            if strcmp(headers{iHeader},'contrast')
                edited_sync_data.contrast = trialArray(:,iHeader);
            elseif strcmp(headers{iHeader},'sound_bit')
                edited_sync_data.sound_bit = trialArray(:,iHeader);
            elseif strcmp(headers{iHeader},'type')
                edited_sync_data.type = trialArray(:,iHeader);
                edited_sync_data.stimIDs = makeUniqueStimID(edited_sync_data.contrast,edited_sync_data.sound_bit,edited_sync_data.type);
            elseif strcmp(headers{iHeader},'sound_dir')
                edited_sync_data.sound_dir = trialArray(:,iHeader);
            end
        end
        
        edited_sync_data.photodiode = photodiode; % doesn't seem necessary to edit this
        
    case 'ori'
        
        % there are four pulses per stimulus
        % edited_sync_data.photodiode = photodiode(1:4:end); % apparently this worked for one experiment
        
        IDs = find(diff(photodiode)>5);
        firstTimes(1) = photodiode(1);
        firstTimes(2:length(IDs)+1) = photodiode(IDs+1);
        edited_sync_data.photodiode = firstTimes';
        edited_sync_data.stimIDs = trialArray(:,2);
        
    case 'retinotopy'
        
        % get stimulus information out of trial array
        stim_x = trialArray(:,2);
        stim_y = trialArray(:,3);
        stim_xy = [stim_x,stim_y];
        
        unique_x = unique(stim_x);
        unique_y = unique(stim_y);
        [p,q] = meshgrid(unique_x, unique_y);
        positions = [p(:) q(:)]; % all possible positions
        
        
        edited_sync_data.stimIDs = stim_xy;
        % reformat sync data
        start_pulse = photodiode(1);
        photodiode_cut(1) = start_pulse;
        
        % this section is stupid but it's just to find the first photodiode
        % pulse for the retinotopy expt
        for i = 1:length(stim_x)
            allIDs = find(photodiode >= start_pulse + i*5);
            we_good = 0;
            for iID = 1:length(allIDs)
                diffID = photodiode(allIDs(iID)) - photodiode_cut(i);
                if we_good == 0
                    if diffID > 3
                        photodiode_cut(i+1) = photodiode(allIDs(iID));
                        we_good = 1;
                    end
                end
            end
            clear allIDs
        end
        
        edited_sync_data.photodiode = photodiode_cut';
end




