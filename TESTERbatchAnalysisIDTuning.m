%load group file!
load ('FullTuningAndLaserAnalysis28-Sep-2016.mat')

fieldNames = fieldnames(fullMaster);
fieldNum = length(fieldNames);

masterIndex = zeros(1,1); %this will keep track of indices of good cells.


IDCells = struct;
IDLatencies = struct;

masterCounter = 1;
%This goes through and 
for i = 1:fieldNum
    numTrodes = fullMaster.(fieldNames{i}).NumberTrodes;
    laserData = fullMaster.(fieldNames{i}).LaserData;
    cellFinder = strfind(fieldnames(laserData),'nt');
    cellFinder = find(not(cellfun('isempty', cellFinder)));
    cellNames = fieldnames(laserData);
    cellNames = cellNames(cellFinder);
    if ~iscell(cellNames)
        cellNames = {cellNames};
    end
    goodCells = cell(1);
    goodCellCount = 1;
    for j = 1:numTrodes
        numClusters = fullMaster.(fieldNames{i}).(cellNames{j}).Clusters;
        goodClusters = [];
        for k = 1:numClusters
            zScoreCross = laserData.(cellNames{j}).zScoreCrossing{k};
            waveCorr = corrcoef(laserData.(cellNames{j}).AverageWaveFormLaser(:,k,2),laserData.(cellNames{j}).AverageWaveFormNonLaser(:,k,2));
            if ~isempty(zScoreCross) & waveCorr(2) > 0.95 & zScoreCross > 0.005
                IDCells.(fieldNames{i}).LaserData.(cellNames{j}) = fullMaster.(fieldNames{i}).LaserData.(cellNames{j});
                IDCells.(fieldNames{i}).(cellNames{j}) = fullMaster.(fieldNames{i}).(cellNames{j});
                infoName = strcat(fieldNames{i},cellNames{j},num2str(k),'StatLatencies');
                IDLatencies.(infoName) = squeeze(fullMaster.(fieldNames{i}).(cellNames{j}).ResponseStatsGraph(k,:,:,:));
                infoName = strcat(fieldNames{i},cellNames{j},num2str(k),'MeanLatencies');
                IDLatencies.(infoName) = squeeze(fullMaster.(fieldNames{i}).(cellNames{j}).FirstSpikeStats(:,:,k,:));
                goodClusters = [goodClusters;k];
                goodCells{goodCellCount} = cellNames{j};
                goodCellCount = goodCellCount + 1;
                masterIndex(masterCounter,1) = i;
                masterIndex(masterCounter,2) = j;
                masterIndex(masterCounter,3) = k;
                masterIndex(masterCounter,4) = zScoreCross;
                masterIndex(masterCounter,5) = waveCorr(2);
                masterCounter = masterCounter + 1;
            end
        end
        if isfield(IDCells,(fieldNames{i}))
            IDCells.(fieldNames{i}).SoundData = fullMaster.(fieldNames{i}).SoundData;
            IDCells.(fieldNames{i}).LaserData.(cellNames{j}).GoodClusters = goodClusters;
            IDCells.(fieldNames{i}).LaserData.(cellNames{j}).GoodCells = goodCells;
            IDCells.(fieldNames{i}).(cellNames{j}).GoodClusters = goodClusters;
        end
    end
    if isfield(IDCells,(fieldNames{i}))
        IDCells.(fieldNames{i}).GoodCells = goodCells;
    end
end

for i = 1:size(masterIndex,1)
    %want to set the path to the correct files for plotting
    
end


save 