function mergeGLXBinFiles(files,nChansTotal,varargin)

%% INPUTS
% files =  {'NP8\20180301\NP8_003_g0_t0.imec.ap','NP8\20180301\NP8_003_2_g0_t0.imec.ap'}
% for opt4: nChansTotal = 277; for opt1,2,3: nChansTotal = 385
% (exclude the '.bin' extension on filenames)


%% OPTIONS
chunkSize = 10000000;

% change the rootDataDir based on where the experiment lives on the server
%rootDataDir = 'B:\ajuavine\Neuropixels\NP9\NP9_003\Raw\'; % grid-hs server
rootDataDir = 'A:\data\ashley_looming\Neuropixels\NP16\NP16_004\'; %nlsas server



if ~isempty(varargin)
    saveFilename = varargin{1};
else
    saveFilename = [rootDataDir files{1} '_merged']; % file will be named after the first
end

%% SET UP FILES
fidOut = [];
fidOut = fopen([saveFilename '.bin'], 'w');

numFiles = length(files);

for iFile = 1:numFiles
    
    % open and save first file
    filename = [rootDataDir files{iFile}];
    d = dir([filename '.bin']);
    nSampsTotal = d.bytes/nChansTotal/2;
    nChunksTotal = ceil(nSampsTotal/chunkSize);
    
    fid = [];
    fid = fopen([filename '.bin'], 'r');
    
    chunkInd = 1;
    
    while 1
        
        fprintf(1, 'file %d/%d - chunk %d/%d\n', iFile, numFiles, chunkInd, nChunksTotal);
        dat = fread(fid, [nChansTotal chunkSize], '*int16');
        
        if ~isempty(dat)
            fwrite(fidOut, dat, 'int16');
        else
            break;
        end
        chunkInd = chunkInd+1;
    end
end

if ~isempty(fid)
    fclose(fid);
end

if ~isempty(fidOut)
    fclose(fidOut);
end

