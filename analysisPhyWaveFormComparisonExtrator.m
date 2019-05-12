
clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);


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

%parameters
masterHeaderSize = 12; %only want the first 12 values of masterData. 
cellCount = 1;
%actually extract files.
for i = 1:numFiles
    %load the target file
    disp(strcat('Analyzing File:',num2str(i)))
    load(targetFiles{i})
    numUnits = size(masterData,1);
    %now lets pull the spikes
    for j = 1:numUnits
        %matclust!
        tarSpikes = squeeze(s.AllSpikesMatclust(:,:,j));
        [C maxInd] = max(max(tarSpikes));
        tarWave = tarSpikes(:,maxInd);
        interpVect = [1:0.1:40];
        interpWave = interp1(1:40,tarWave,interpVect,'spline');
        
        %get peak trough
        [pkVal pkInd] = max(interpWave(100:end));
        pkInd = pkInd + 100 - 1;
        %now we need to find the minimum following the peak
        [troughVal troughInd] = min(interpWave(pkInd:end));
        troughInd = troughInd + pkInd - 1;
        peakTrough(cellCount,1) = (troughInd - pkInd)/300000;
        
        %get half-width. Find before peak. 
        firstFind = find(interpWave(1:pkInd) > pkVal/2,1,'first');
        %now get the second. 
        secondFind = find(interpWave(firstFind:end) < pkVal/2,1,'first');
        halfWidth(cellCount,1) = secondFind/300000;
        
        waveStoreMat(cellCount,:) = interpWave;
        
        %now pull from phy
        tarSpikes = squeeze(s.AllSpikesPhy(j,:,:));
        tarSpikes = tarSpikes - tarSpikes(:,1);
        %we need to eliminate big drifties. To do so, check last value. 
        findBigNeg = find(tarSpikes(:,end) < -500);
        tarSpikes(findBigNeg,:) = 0;
        [c minInd] = min(squeeze(tarSpikes(:,20)));
        tarWave = tarSpikes(minInd,:);
        %apply mild smoothing
        tarWave = smooth(tarWave,3);
        %interpolate
        interpVect = [1:0.1:54];
        interpWave = interp1(1:54,tarWave,interpVect,'spline');
        
        %get peak trough, start around equivalent of 20
        [pkVal pkInd] = min(interpWave(100:end));
        pkInd = pkInd + 100 - 1;
        %now we need to find the minimum following the peak
        [troughVal troughInd] = max(interpWave(pkInd:end));
        troughInd = troughInd + pkInd - 1;
        peakTrough(cellCount,2) = (troughInd - pkInd)/300000;
        
        %get half-width. Find before peak. 
        firstFind = find(interpWave(1:pkInd) < pkVal/2,1,'first');
        %now get the second. 
        secondFind = find(interpWave(firstFind:end) > pkVal/2,1,'first');
        halfWidth(cellCount,2) = secondFind/300000;
        
        waveStorePhy(cellCount,:) = interpWave;
        
        cellCount = cellCount + 1;
        
    end
    
end


figure
hold on
plot(peakTrough(:,1),halfWidth(:,1),'k.')
plot(peakTrough(:,2),halfWidth(:,2),'r.')
for i = 1:length(peakTrough)
    quiver(peakTrough(i,1),halfWidth(i,1),peakTrough(i,2)-peakTrough(i,1),halfWidth(i,2)-halfWidth(i,1),0)
end


figure
hold on
plot(peakTrough(:,1),halfWidth(:,1),'ko')
plot(peakTrough(:,2),halfWidth(:,2),'ro')


findPVs = find(peakTrough(:,2) < 0.00055);












