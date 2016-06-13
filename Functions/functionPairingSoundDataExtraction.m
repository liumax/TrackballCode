%This function is meant to calculate information of tuning curves based on
%selected subset of spikes. 

%NOTE: spikeName and ttlName must be in quotations.
function [masterStruct] = functionPairingSoundDataExtraction(masterStruct,...
    soundName,soundFile); 
%pulls out appropriate sound file.
soundData = soundFile.(soundName);
%opens empty array for storage of values
master = zeros(size(soundData.Frequencies,1),5);
%saves all frequencies delivered, figures out unique frequencies delivered
master(:,2) = soundData.Frequencies;
uniqueFreqs = unique(master(:,2));
%saves all dBs delivered, figures out unique dBs delivered
master(:,3) = soundData.dBs;
uniqueDBs = unique(master(:,3));
%generates easy placeholders for later for loops
numFreqs = size(uniqueFreqs,1);
numDBs = size(uniqueDBs,1);

%% Finds the number of octaves, makes array of octave steps. This will be used for imagesc graphing applications
totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
%This then makes an array of the full octave steps I've made
octaveRange = zeros(totalOctaves + 1,2);
octaveRange(1,1) = uniqueFreqs(1);
for i = 1:totalOctaves
    octaveRange (i+1,1) = octaveRange(i,1)*2;
end
%next, I find the positions from uniqueFreqs that match octaveRange
for i = 1:size(octaveRange,1);
    octaveRange(i,2) = find(round(uniqueFreqs) == octaveRange(i,1));
end

%% Does the same for dBs. 
dbSteps = uniqueDBs(2) - uniqueDBs(1);
totalDBs = (uniqueDBs(end) - uniqueDBs(1))/dbSteps;
dbRange = zeros(totalDBs + 1,2);
dbRange(:,1) = uniqueDBs(1):dbSteps:uniqueDBs(end);
for i = 1:size(dbRange,1)
    dbRange(i,2) = find(uniqueDBs == dbRange(i,1));
end

%% Generates index that goes from low freq to high freq, and within each freq, goes from low to high amplitude
master(:,4) = 1:1:size(master,1);
master(:,5) = zeros;

sortingCounter = 1;

for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == uniqueDBs(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end
%now, master(:,5) is the index if I want to sort rasters by frequency and
%amplitude

%% Generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end

%% Save information into masterStructure
masterStruct.SoundData.(soundName).UniqueFreqs = uniqueFreqs;
masterStruct.SoundData.(soundName).UniqueDBs = uniqueDBs;
masterStruct.SoundData.(soundName).Frequencies = master(:,2);
masterStruct.SoundData.(soundName).dBs = master(:,3);
masterStruct.SoundData.(soundName).ToneReps = soundData.ToneRepetitions;
masterStruct.SoundData.(soundName).ToneDur = soundData.ToneDuration;
masterStruct.SoundData.(soundName).MasterArray = master;
masterStruct.SoundData.(soundName).OctaveRange = octaveRange;
masterStruct.SoundData.(soundName).dBRange = dbRange;
masterStruct.SoundData.(soundName).FrequencyNames = freqNameHolder;

end