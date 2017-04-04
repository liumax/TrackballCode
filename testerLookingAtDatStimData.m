











%this is meant to be master code to run analysis code!
%Establishes the home folder as the base
masterFolder = pwd;
masterDir = dir;
masterDir = {masterDir.name};
%this code here is to remove folders reported in dir like "." or ".."
masterIndex = strfind(masterDir,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
%masterFolders is a listing of all folders that I want to target for batch
%analysis
masterFolders = masterDir(masterIndex);
numFolders = size(masterFolders,2);
for masterCount = 1:numFolders
%generates text that signifies path to the target folder
targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
cd(targetPath)
%what lists matlab files in the folder. can extract based on file type.
test = what;
test = test.mat;
test = test{1};
load(test);
nameFind = strfind(test,'DatStimAnalysis');
newName =strcat('x',test(1:nameFind-1));
fullData.(newName) = s;
cd(masterFolder)
end


%pull names of recordings out
fieldNames = fields(fullData);
%pull individual components of names out
fieldData = cell(length(fieldNames),5);
for i = 1:length(fieldNames)
    %pull name as string
    nameHold = fieldNames{i};
    spaceFind = strfind(nameHold,'_');
    %store date
    fieldData{i,1} = nameHold(2:spaceFind(1)-1);
    %store animal name
    fieldData{i,2} = nameHold(spaceFind(1)+1:spaceFind(2)-1);
    %store hemisphere
    fieldData{i,3} = nameHold(spaceFind(2)+1:spaceFind(3)-1);
    %store depth
    fieldData{i,4} = nameHold(spaceFind(3)+1:spaceFind(4)-1);
    %store type
    fieldData{i,5} = nameHold(spaceFind(4)+1:end);
end




store.FileName = cell(100,1);
store.UnitName = cell(100,1);
store.BaselineRate = zeros(100,1);
store.BaselineSTD = zeros(100,1);
store.PostZeroBin = zeros(100,1);
store.PostZeroSTD = zeros(100,1);
store.PostBigBin = zeros(100,1);
store.TrueAUC = zeros(100,1);
store.ShuffleAUC = zeros(100,2);
store.Waves = zeros(40,4,100);
placeHolder = 1;
for i = 1:length(fieldNames)
    %load data
    testData = fullData.(fieldNames{i});
    desigNames = testData.DesignationName;
    %pull tuning data
    %go into individual units to pull data
    for j = 1:length(desigNames)
        store.FileName{placeHolder - 1 + j} = fieldNames{i};
        store.UnitName{placeHolder - 1 + j} = desigNames{j};
        store.StartPref(placeHolder - 1 + j) = testData.RotaryData.StartPreference;
        store.EndPref(placeHolder - 1 + j) = testData.RotaryData.EndPreference;
        store.PrefRange((placeHolder - 1 + j),1) = prctile(testData.RotaryData.PreferenceDistribution,0.5);
        store.PrefRange((placeHolder - 1 + j),2) = prctile(testData.RotaryData.PreferenceDistribution,99.5);
        store.BaselineRate(placeHolder - 1 + j) = testData.(desigNames{j}).AverageRate;
        store.BaselineSTD(placeHolder - 1 + j) = testData.(desigNames{j}).AverageSTD;
        store.PostZeroBin(placeHolder - 1 + j) = testData.(desigNames{j}).AllHistograms(81);
        store.PostZeroSTD(placeHolder - 1 + j) = testData.(desigNames{j}).HistogramStandardDeviation(81);
        store.PostBigBin(placeHolder - 1 + j) = mean(testData.(desigNames{j}).AllHistograms(81:101));
        if isfield(testData.(desigNames{j}),'TrueAUC')
        store.TrueAUC(placeHolder - 1 + j) = testData.(desigNames{j}).TrueAUC;
        store.ShuffleAUC(placeHolder - 1 + j,1) = prctile(testData.(desigNames{j}).ShuffleAUC,0.5);
        store.ShuffleAUC(placeHolder - 1 + j,2) = prctile(testData.(desigNames{j}).ShuffleAUC,99.5);
        end
        %store waveform
        store.Waves(:,:,placeHolder - 1 + j) = testData.(desigNames{j}).AverageWaveForms;
    end
    placeHolder = placeHolder + length(desigNames);
end

%lets plot some basics

figure
hist(store.BaselineRate,100);
title('Distribution of Baseline Firing Rates')

for i = 1:placeHolder - 1
    if store.TrueAUC(i) > store.ShuffleAUC(i,2) | store.TrueAUC(i) < store.ShuffleAUC(i,1)
        findSigAUC(i) = 1;
    else
        findSigAUC(i) = 0;
    end
end
findSigAUC = find(findSigAUC == 1);

figure
plot(store.TrueAUC,store.BaselineRate,'k.')
hold on
plot(store.TrueAUC(findSigAUC),store.BaselineRate(findSigAUC),'ro')
title({'AUC Values vs Baseline Rate';'Red Circles are Significant by 99% CI'})


%okay, now lets get at the running. DO we see preference for running or
%stopping?

for i = 1:placeHolder - 1
    if store.StartPref(i) > store.PrefRange(i,2) | store.StartPref(i) < store.PrefRange(i,1)
        findSigStart(i) = 1;
    else
        findSigStart(i) = 0;
    end
end
findSigStart = find(findSigStart == 1);

















