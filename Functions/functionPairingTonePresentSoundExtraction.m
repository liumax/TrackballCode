%This function is meant to extract sound information from tone
%presentations without pairing.

%NOTE: spikeName and ttlName must be in quotations.
function [masterStruct] = functionPairingTonePresentSoundExtraction(masterStruct,...
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

%% Generates index that goes from low freq to high freq, and within each freq, goes from low to high amplitude
master(:,4) = 1:1:size(master,1);

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
masterStruct.SoundData.(soundName).FrequencyNames = freqNameHolder;
masterStruct.SoundData.(soundName).TargetSet = [soundData.TargetFreq,soundData.TargetDB];
masterStruct.SoundData.(soundName).ControlSet = [soundData.ControlFreq,soundData.ControlDB];

end