%load group file!
load ('FullTuningAndLaserAnalysis28-Sep-2016.mat')

fieldNames = fieldnames(fullMaster);
fieldNum = length(fieldNames);

% masterIndex = zeros(1,1); %this will keep track of indices of good cells.
masterInfo = zeros(1,1);
histHolder = cell(1,1);
masterCounter = 1;
%This goes through and 
for i = 1:fieldNum
    numTrodes = fullMaster.(fieldNames{i}).NumberTrodes;
    cellFinder = strfind(fieldnames(fullMaster.(fieldNames{i})),'nt');
    cellFinder = find(not(cellfun('isempty', cellFinder)));
    cellNames = fieldnames(fullMaster.(fieldNames{i}));
    cellNames = cellNames(cellFinder);
    if ~iscell(cellNames)
        cellNames = {cellNames};
    end
    for j = 1:numTrodes
        RPVs= fullMaster.(fieldNames{i}).(cellNames{j}).RPVs;
        goodRPVs = find(RPVs < 0.1);
        for k = 1:size(goodRPVs,1)
            masterInfo(masterCounter,1) = fullMaster.(fieldNames{i}).(cellNames{j}).AverageRate(k);
            histHolder{masterCounter,1} = squeeze(fullMaster.(fieldNames{i}).(cellNames{j}).AllHistograms(:,k));
            histHolder{masterCounter,2} = squeeze(fullMaster.(fieldNames{i}).(cellNames{j}).ResponseStatsGraph(k,:,:,:));
            masterCounter = masterCounter + 1;
            plot(squeeze(fullMaster.(fieldNames{i}).(cellNames{j}).AllHistograms(:,k)))
        end
    end
end



for i = 1:size(masterIndex,1)
    %want to set the path to the correct files for plotting
    
end


save 