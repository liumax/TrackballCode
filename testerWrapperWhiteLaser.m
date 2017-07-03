% testerWrapperWhiteLaser

%pull the analysis files that have already been processed.
targetFiles = {'170629_ML170621G_R10_2800_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat'};

%get the number of files
numFiles = length(targetFiles);

%extract distances
for i = 1:numFiles
    targetName = targetFiles{i};
    findSpace = strfind(targetName,'_');
    probeDepth(i) = str2num(targetName(findSpace(3)+1:findSpace(4)-1));
end

%now we need to make a big master sheet.
masterInd = 1;
for i = 1:numFiles
    %load the target file
    load(targetFiles{i})
    %extract master array
    master = s.MasterSheet;
    masterHeader = s.MasterDesigs;
    %pull indices
    [indPVMSN] = functionCellStringFind(masterHeader,'PV/MSN Desig');
    [indPkTrough] = functionCellStringFind(masterHeader,'PeakTroughTimeinMS');
    [indOverFire] = functionCellStringFind(masterHeader,'OverallFiringRate');
    [indDistance] = functionCellStringFind(masterHeader,'Distance from Top of Shank');
    [indShankDes] = functionCellStringFind(masterHeader,'Shank Designation');
    [indPreAverage] = functionCellStringFind(masterHeader,'PreAverage');
    [indPostAverage] = functionCellStringFind(masterHeader,'PostAverage');
    [indLaserOnly] = functionCellStringFind(masterHeader,'LaserOnlyAverage');
    [indToneOnly] = functionCellStringFind(masterHeader,'ToneOnlyAverage');
    [indToneLaser] = functionCellStringFind(masterHeader,'ToneLaserAverage');
    [indPreAverageRunning] = functionCellStringFind(masterHeader,'Pre AverageRunning');
    [indPreAverageStationary] = functionCellStringFind(masterHeader,'Pre AverageStationary');
    [indLocoAUC] = functionCellStringFind(masterHeader,'LocoAUCScore');
    [indLocoSig] = functionCellStringFind(masterHeader,'LocoAUCSignificance');
    
    %now i need to adjust the depth system. First, add the height of the
    %probe in terms of number of probe sites. then multiply by negative 1 and 12
    master(:,indDistance) = (master(:,indDistance) + s.ShankLength)*-12.5 + probeDepth;
    numUnits = size(master,1);
    fullMaster(masterInd: masterInd + numUnits - 1,:) = master;
    fullMasterHeaders(:,i) = masterHeader;
    masterInd = masterInd + numUnits;
    
    
end













