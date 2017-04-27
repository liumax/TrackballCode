%This is meant to be code that will take the data from my combined dataset
%and process it into something more amenable for analysis and comparisons
%between groups


%pull names of recordings out
fieldNames = fields(fullMaster);
%pull individual components of names out
fieldData = cell(length(fieldNames),5);
counter = 1;

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
    unitPer(i) = length(fullMaster.(fieldNames{i}).DesignationName);
    recordNames(counter:counter+unitPer(i)-1,1) = i;
    counter = counter + unitPer(i);
end

%now we need to go through and sort by animal and hemisphere. Each
%hemisphere should be selectable on its own.
hemiDesig = zeros(length(fieldNames),1);
counter = 1;

for i = 1:length(fieldNames)-1
    %pull file and hemisphere designations
    currFile = fieldData{i,2};
    currHemi = fieldData{i,3};
    nextFile = fieldData{i+1,2};
    nextHemi = fieldData{i+1,3};
    %look for match in file name
    if strcmp(currFile,nextFile)
        if strcmp(currHemi,nextHemi)
            hemiDesig(i) = counter;
            hemiDesig(i+1) = counter;
        else
            hemiDesig(i) = counter;
            counter = counter + 1;
            hemiDesig(i+1) = counter;
        end
    else
        hemiDesig(i) = counter;
        counter = counter + 1;
        hemiDesig(i+1) = counter;
    end
end

%put the experiment number into recordNames. 
recordNames(:,2) = zeros(length(recordNames),1);
for i = 1:length(unique(hemiDesig))
    %find the file names associated with the target hemisphere
    hemiFind = find(hemiDesig == i);
    recordNames(ismember(recordNames,hemiFind),2) = i;
end
%make big structured array for storage of recording info
for i = 1:length(recordNames)
    infoStruct(i).Date = fieldData{recordNames(i,1),1};
    infoStruct(i).Animal = fieldData{recordNames(i,1),2};
    infoStruct(i).Hemisphere = fieldData{recordNames(i,1),3};
    infoStruct(i).Depth = fieldData{recordNames(i,1),4};
    infoStruct(i).Condition = fieldData{recordNames(i,1),5};
end

%now lets pull PV vs MSN identity, as well as AUC scores and significance
%while performing other stashing. 
unitStore = struct;
recordNames(:,3:4) = zeros(length(recordNames),2);
counter = 1;
for i = 1:length(fieldNames)
    tempData = fullMaster.(fieldNames{i});
    desigNames = tempData.DesignationName;
    numUnits = length(desigNames);
    for j = 1:numUnits
        %fill in information structure
        infoStruct(counter).UnitID = desigNames{j};
        %fill in the structured table
        unitStore(counter).UnitID = desigNames{j};
        unitStore(counter).TrialNumber = size(tempData.(desigNames{j}).IndividualHistograms,2);
        unitStore(counter).AverageWaves = tempData.(desigNames{j}).AverageWaveForms;
        unitStore(counter).AverageHistogram = hist(tempData.(desigNames{j}).AllRasters(:,1),tempData.(desigNames{j}).HistBinVector)/unitStore(counter).TrialNumber/tempData.Parameters.histBin;
        unitStore(counter).HistogramVector = tempData.(desigNames{j}).HistBinVector;
        unitStore(counter).Rasters = tempData.(desigNames{j}).AllRasters;
        unitStore(counter).RasterWindow = tempData.Parameters.RasterWindow;
        unitStore(counter).AllSpikes = tempData.(desigNames{j}).SpikeTimes;
        unitStore(counter).AUCScore = tempData.(desigNames{j}).TrueAUC;
        unitStore(counter).AUCRange = [prctile(tempData.(desigNames{j}).ShuffleAUC,0.05),prctile(tempData.(desigNames{j}).ShuffleAUC,9.95)];
        
        
        %lets pick out laser times!
        laserOn = tempData.LaserData.LaserStartTimes;
        laserOff = tempData.LaserData.LaserEndTimes;
        laserDur = tempData.LaserData.LaserDuration;
        %figure out pulses per trial
        pulsePerTrial = round(length(laserOn)/unitStore(counter).TrialNumber);
        unitStore(counter).LaserStarts = laserOn(1:pulsePerTrial:end);
        unitStore(counter).LaserEnds = laserOff(pulsePerTrial:pulsePerTrial:end);
        %store laser duration and pulse information
        infoStruct(counter).TrialNumber = unitStore(counter).TrialNumber;
        infoStruct(counter).LaserPulses = pulsePerTrial;
        infoStruct(counter).LaserDuration = laserDur;
        aveLaserDur = mean(unitStore(counter).LaserEnds - unitStore(counter).LaserStarts);
        
        %lets do some calculations!
        %first, calculate baseline firing rate using the baseline period
        %of rasterWindow(1) : 0 from raster data
        for k = 1: unitStore(counter).TrialNumber
            %see if any spikes in this trial
            spikeFinder = length(find(unitStore(counter).Rasters(:,2) == k & unitStore(counter).Rasters(:,1) <0));
            %calculate baseline rate
            baseRate(k) = spikeFinder/(-1*unitStore(counter).RasterWindow(1));
            %now find spikes in the laser period
            spikeFinder = length(find(unitStore(counter).Rasters(:,2) == k & unitStore(counter).Rasters(:,1) <= aveLaserDur & unitStore(counter).Rasters(:,1) > 0));
            laserRate(k) = spikeFinder/(aveLaserDur);
            %now find spikes in the postlaser period
            spikeFinder = length(find(unitStore(counter).Rasters(:,2) == k & unitStore(counter).Rasters(:,1) <= aveLaserDur + 3 & unitStore(counter).Rasters(:,1) > aveLaserDur));
            postLaserRate(k) = spikeFinder/(3);
        end
        unitStore(counter).BaselineRate = mean(baseRate);
        unitStore(counter).BaselineSTD = std(baseRate);
        %now lets make a z-scored histogram
        unitStore(counter).zHist = (unitStore(counter).AverageHistogram - unitStore(counter).BaselineRate)/unitStore(counter).BaselineSTD;
        %lets also make z-scores of different periods
        unitStore(counter).zLaserMean = (mean(laserRate)- unitStore(counter).BaselineRate)/unitStore(counter).BaselineSTD;
        unitStore(counter).zPostLaserMean = (mean(postLaserRate)- unitStore(counter).BaselineRate)/unitStore(counter).BaselineSTD;
        %also store raw values
        unitStore(counter).BaselineCounts = baseRate;
        unitStore(counter).LaserCounts = laserRate;
        unitStore(counter).PostLaserCounts = postLaserRate;
        clear 'baseRate' 'laserRate' 'postLaserRate'
        
        %lets do some analysis of waveform to classify
        targetWaves = unitStore(counter).AverageWaves;
        [maxVal maxFind] = max(max(targetWaves));
        sortedPeakValues = maxVal;
        sortedWaves = interp1([1:1:40],targetWaves(:,maxFind),[1:0.1:40],'spline');
        %turns out this seems to work fairly well. now lets pull peaks and
        %troughs
        %the peak should be the biggest point in the whole thing. Lets pull
        %that out
        [mVal,sortedPeaks] = max(sortedWaves(100:200));
        %the trough should be the opposite
        [minVal,sortedPeakTroughs] = min(sortedWaves(sortedPeaks:end));
        sortedPeakTroughs = sortedPeakTroughs/300000;
        unitStore(counter).PeakTrough = sortedPeakTroughs;
        unitStore(counter).OverallRate = tempData.(desigNames{j}).OverallFiringRate;
        %classify based on this value
        if sortedPeakTroughs > 0.0008
            recordNames(counter,3) = 1;
            infoStruct(counter).CellClass = 'MSN';
        elseif sortedPeakTroughs < 0.0008 && tempData.(desigNames{j}).OverallFiringRate > 2
            recordNames(counter,3) = 2;
            infoStruct(counter).CellClass = 'PV';
        else
            infoStruct(counter).CellClass = 'UNK';
        end
        
        %now lets do some basic AUC analysis
        if unitStore(counter).AUCScore > unitStore(counter).AUCRange(2)
            infoStruct(counter).AUCClass = '+';
            recordNames(counter,4) = 1;
        elseif unitStore(counter).AUCScore < unitStore(counter).AUCRange(1)
            infoStruct(counter).AUCClass = '-';
            recordNames(counter,4) = -1;
        else
            infoStruct(counter).AUCClass = '0';
        end
        counter = counter + 1;
    end
end

recordKey = {'Index of File Name','Index of Hemispheric Recording','Cell ID (1 = MSN, 2 = PV)','AUC Significance'};

save('dataset.mat','infoStruct','recordNames','unitStore','recordKey')