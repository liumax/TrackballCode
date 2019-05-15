%% This is solely to extract units suitable for MID analysis from the new dataset. (after 4k switch)

clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'Analysis');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);

minSpikes = 1000;
bigStoreSpikes = [];
bigStoreWaves = [];
bigStorePkTrough = [];
bigStoreSTA= [];
bigStoreCellType = [];
masterInd = 1;

%parameters
masterHeaderSize = 12; %only want the first 12 values of masterData. 

%actually extract files.
for i = 1:numFiles
    %load the target file
    disp(strcat('Analyzing File:',num2str(i)))
    load(targetFiles{i})
    numCells = length(s.DesignationName);
    desigName = s.DesignationName;
    
    %now lets pull significance information from STA
    newTarget = targetFiles{i};
    newTarget = newTarget(1:end-22);
    newTarget = [newTarget,'DMRData.mat'];
    load(newTarget)
    
    realCorrStore = realCorrStore(1:numCells);
    corrCoefStore = corrCoefStore(1:numCells,:);
    staSpikes = sum(spikeArray');
    
    for j = 1:numCells
        disp(desigName{j})
        %determine significance of STA. 
        corrPrctile = prctile(corrCoefStore(j,:),[0.5 99.5]);
        %determine number of spikes. 
        if realCorrStore(j) > corrPrctile(2) & staSpikes(j) > minSpikes %significant STA with enough spikes
            disp('Significant STA Found')
            disp([targetFiles{i},desigName{j},])
            %store name
            nameStore{masterInd} = targetFiles{i};
            unitStore{masterInd} = desigName{j};
            %store cell type
            bigStoreCellType(masterInd) = masterData(j,7);
            %store spikes
            bigStoreSpikes(masterInd,:) = spikeArray(j,:);
            %store STA
            bigStoreSTA(masterInd,:) = s.STAs(j,:);
            %pull up cell average waveforms
            cellWaves = s.(desigName{j}).AverageWaveForms;
            %now pull up which one is biggest
            waveMax = max(cellWaves);
            [maxVal maxInd] = max(waveMax);
            %now generate finely sampled version
            chosenWave = cellWaves(:,maxInd);

            interpVect = [1:0.1:40];
            interpWave = interp1(1:40,chosenWave,interpVect,'spline');
            %Store wave
            bigStoreWaves(masterInd,:) = interpWave;
            %now find peak value
            [pkVal pkInd] = max(interpWave(100:end));
            pkInd = pkInd + 100 - 1;
            %now we need to find the minimum following the peak
            waveDiff = diff(interpWave);

            troughInd = find(waveDiff(pkInd:end) > 0,1,'first');
            troughVal = interpWave(troughInd);
            if isempty(troughInd)
                [troughVal troughInd] = min(interpWave(pkInd:end));
            end
            troughInd = troughInd + pkInd - 1;
            
            bigStorePkTrough(masterInd) = pkVal/troughVal;
            
            
            
            
            masterInd = masterInd + 1;
            
        end
    end
end