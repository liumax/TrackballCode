%This is code that is meant to run through a bunch of folders and extract
%average LFP traces for each nTrode aligned to DIO1. 

%pull out directory names
dirCont = dir;
dirNames = extractfield(dirCont,'name');
isDirList = cell2mat(extractfield(dirCont,'isdir'));
trueDirs = find(isDirList);

homeDir = pwd;

lfpStore = struct;
rasterWindow = [-4,6];

counter = 1;
nameArray = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W'};
for i = 3:length(trueDirs); %use 3 to avoid the first two fake folders
    disp((dirNames{trueDirs(i)}))
    %get into targeted folder
    cd (dirNames{trueDirs(i)})
    %pull DIO times
    x = dir;
    xNames = extractfield(x,'name');
    folderFinder = strfind(xNames,'DIO');
    folderFinder = find(~cellfun(@isempty,folderFinder));
    targetFolder = xNames{folderFinder};
    %enter target folder
    cd (targetFolder)
    %find target file
    x = dir;
    xNames = extractfield(x,'name');
    folderFinder = strfind(xNames,'D1.dat');
    folderFinder = find(~cellfun(@isempty,folderFinder));
    targetFile = xNames{folderFinder};
    %extract data
    [DIOData] = readTrodesExtractedDataFile(targetFile);
    [dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,30000);
    %for dat stim, dioTimes are my targets!
    
    %now return to main folder
    cd ..
    %and find LFP folder
    x = dir;
    xNames = extractfield(x,'name');
    folderFinder = strfind(xNames,'LFP');
    folderFinder = find(~cellfun(@isempty,folderFinder));
    targetFolder = xNames{folderFinder};
    %enter LFP folder
    cd (targetFolder)
    %pull individual file names
    x = dir;
    xNames = extractfield(x,'name');
    folderFinder = strfind(xNames,'LFP');
    folderFinder = find(~cellfun(@isempty,folderFinder));
    targetNames = xNames(folderFinder);
    %now go into individual LFPs
    for j = 1:length(targetNames)
        lfp = readTrodesExtractedDataFile(targetNames{j});
        viewSamples = round((rasterWindow(2)-rasterWindow(1))*lfp.clockrate/lfp.decimation);
        %counts number of LFP samples
        lfpSamples = size(lfp.fields.data,1);
        %makes LFP time points based on # of samples and decimation
        %adjust time to actual time (to match master)
        lfpTimes = (((0:1:lfpSamples-1)*lfp.decimation)+lfp.first_timestamp)'/30000;
        lfpSignals = lfp.fields.data;
        lfpHolder = [];
        for k = 1:length(dioTimes)
            finder = find(lfpTimes > dioTimes(k) + rasterWindow(1),1);
            lfpHolder(:,k) = lfpSignals(finder:finder+viewSamples-1);
        end
        lfpStore.(nameArray{counter})(j).RawLFP = lfpHolder([1:10:end],:);
        lfpStore.(nameArray{counter})(j).AverageLFP = mean(lfpStore.(nameArray{counter})(j).RawLFP,2);
        lfpStore.(nameArray{counter})(j).StdLFP = std(lfpStore.(nameArray{counter})(j).RawLFP,0,2);
         
    end
    counter = counter + 1;
    cd ..
    cd ..
end