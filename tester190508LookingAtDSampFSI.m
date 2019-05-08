
clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'Analysis');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);

histLims = [-0.2 0.4];

bigMaster = [];
bigMasterInd = 1;
respVect = [];
pkTroughRatioStore = [];
binValBigStore = [];
sigValBigStore = [];
latMapBigStore = [];
widthLatStore = [];
nameStore = [];
unitStore = [];
halfWidthTimeStore = [];
pkTroughTimeStore = [];
dmrSpikeNum = [];
bigpspk = [];
bigpx = [];
bigpxspk = [];


%parameters
masterHeaderSize = 12; %only want the first 12 values of masterData. 

load('/Users/maxliu/Desktop/dmrDatasets/tempVect.mat')
load('/Users/maxliu/Desktop/dmrDatasets/ttlOnsetTime.mat')

masterCount = 1;
%actually extract files.
for i = 1:numFiles
    %load the target file
    disp(strcat('Analyzing File:',num2str(i)))
    load(targetFiles{i})
    desigNames = s.DesignationName;
    numUnits = size(masterData,1);
    %find FSIs
    cellTypes = masterData(:,7);
    findFSIs = find(cellTypes == 1);
    %now pull information to make hist vector
    dmrDIO = s.SoundData.DMRPulses;
    trueDMRtimes = interp1(ttlOnsetTime,dmrDIO,tempVect);
    dmrStep = mean(mode(diff(trueDMRtimes)));
    newWin = round(s.Parameters.STRFWin/dmrStep);
    spikeHistVect = trueDMRtimes+dmrStep/2;
    
    for j = 1:length(findFSIs)
        %pull spikes
        spikeStore = s.(desigNames{findFSIs(j)}).SpikeTimes(s.(desigNames{findFSIs(j)}).SpikeTimes > spikeHistVect(1) + s.Parameters.STRFWin(1) & s.(desigNames{findFSIs(j)}).SpikeTimes < spikeHistVect(end));
        bigSpikeStore{masterCount} = spikeStore - dmrDIO(1) + ttlOnsetTime(1);
        masterCount = masterCount + 1;
    end
    
end


%and now generate a spike matrix!
numUnits = size(bigSpikeStore,2);
normMatrix = [];
for i = 1:numUnits
    spikeStore = bigSpikeStore{i};
    normMatrix(i,:) = hist(spikeStore,tempVect);
end

%now generate STAs
numLags = 100;
load('/Users/maxliu/Desktop/dmrDatasets/stim4k.mat')
[sta, stabigmat, spkcountvec] = quick_calc_sta(stimulus, normMatrix, numLags);



trueDMRtimes = interp1(ttlOnsetTime,dmrDIO,tempVect);


