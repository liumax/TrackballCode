



%identify files, pull names, set up for loop for extraction

clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);

%generate empty things here/store variables

%variables
sigThresh = 5; %minimum number of significant responses at 0.05 p threshold. 
sigWin = 2; %minimum width for entry to analysis. 
tarWin = 2; %target window of analysis. 1 = fast, 2 = tone, 3 = gen. 
boundsLatHist = [-0.6 0.4];

%counter of number of cells I'm keeping
fullCounter = 1;
fullNames = []; %for storage of file and unit names
fullMaster = [];
fullWidth = []; %for storage of width values
fullBin = [];
fullSig = [];
fullLat = [];
latFineHist = [];

bigMasterInd = 1;

for i = 1:numFiles
    %load the target file
    load(targetFiles{i})
    disp(strcat('Loading File-',num2str(i)))
    
    tempBinStore = [];
    dbStore = [];
    freqStore = [];
    
    numUnits = size(masterData,1);
    desigName = s.DesignationName;
    bigMaster(bigMasterInd:bigMasterInd + numUnits - 1,:) = masterData(:,1:12);
    
    %save frequency values
    freqValStore{i} = s.SoundData.UniqueFrequencies(1+s.SoundData.WhiteNoise:end);
    
    %pull trial matrix, replace true DB values with fixed DB values
    trialMatrix = [];
    trialMatrix = s.SoundData.TrialMatrix;
    
    %now we need to go through and replace the actual real DB values with
    %standardized values.
    for j = 1:s.SoundData.NumFreqs
        finder = find(trialMatrix(:,2) == s.SoundData.UniqueFrequencies(j));
        tarDBs = (trialMatrix(finder,3));
        uniqueTars = unique(tarDBs);
        for k = 1:length(uniqueTars)
            trialMatrix(finder(tarDBs == uniqueTars(k)),3) = s.SoundData.UniqueDBs(k);
        end
    end
    %generate a finder for non-white noise tones
    toneFinder = find(s.SoundData.TrialMatrix(:,2) ~= 0);
    
    %pull trial matrix for laser trials, replace true DB values with fixed DB values
    trialMatrixLaser = [];
    trialMatrixLaser = s.SoundData.TrialMatrixPaired;
    %now we need to go through and replace the actual real DB values with
    %standardized values.
    for j = 1:s.SoundData.NumFreqs
        finder = find(trialMatrixLaser(:,2) == s.SoundData.UniqueFrequencies(j));
        tarDBs = (trialMatrixLaser(finder,3));
        uniqueTars = unique(tarDBs);
        for k = 1:length(uniqueTars)
            trialMatrixLaser(finder(tarDBs == uniqueTars(k)),3) = s.SoundData.UniqueDBs(k);
        end
    end
    %generate a finder for non-white noise tones
    toneFinderLaser = find(s.SoundData.TrialMatrixPaired(:,2) ~= 0);
    
    %now we need to clean up tone time data. 
    gapFind = find(s.SoundData.ToneTimeDiff < 0.1);
    gapFindPlus = gapFind + 1;
    fullGap = [gapFind;gapFindPlus];
    %generate new tone thing. 
    newSound = s.SoundData.ToneTimes;
    newSound(fullGap) = [];
    newSoundNon = newSound(2:2:end);
    newSoundLaser = newSound(1:2:end);
    
    uniqueDBs = s.SoundData.UniqueDBs;
    
    %now go through individual units to extract information
    for j = 1:numUnits
        disp(desigName{j})
        
        
        
        %store db data and freq data
        trueDBs{bigMasterInd + j - 1} = uniqueDBs;
        numDBUnit(bigMasterInd + j - 1) = length(uniqueDBs);
        numFreqUnit(bigMasterInd + j - 1) = length(s.SoundData.UniqueFrequencies)-s.SoundData.WhiteNoise;
        
        %find max response, use this to select flanking frequencies to be
        %used for latency analysis.
        [row,col] = find(s.(desigName{j}).BinTone(1+s.SoundData.WhiteNoise:end,:) == max(max(s.(desigName{j}).BinTone(1+s.SoundData.WhiteNoise:end,:))));
        if length(row) == 1
            disp('Found single max value!')
        elseif length(row) > 1
            disp('Found Multiple Max Values. Pruning to max amplitude')
            col = col(end);
            row = row(end);
        else
            disp('No Max Value Found')
        end

        %find appropriate histogram, store.
        if length(row) == 1
            dbVal = col;
            freqVal= row;
%             tempHist(:,j) =
%             squeeze(s.(desigName{j}).FreqDBHistograms(freqVal+s.SoundData.WhiteNoise,end-5+dbVal,:));
            dbStore(j) = dbVal;
            freqStore(j) = freqVal;
        else
%             tempHist(:,j) =
%             zeros(length(s.(desigName{j}).FreqDBHistograms(1,1,:)),1);
        end
        %now we need to pull the peak values from the 70DB band.
        tarVals = s.(desigName{j}).BinTone(1+s.SoundData.WhiteNoise:end,end);
        ind1 = s.SoundData.NumDBs - 1;
        ind2 = s.SoundData.NumDBs;
        [pk ind] = max(tarVals);
        if ind >=2 & ind <= length(s.SoundData.UniqueFrequencies) - s.SoundData.WhiteNoise - 1 %these are non-edge cases
            tarCols = [ind1,ind1,ind1,ind2,ind2,ind2];
            tarRows = [ind-1+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise,ind-1+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise];
        elseif ind < 2 %edge cases. 
            tarCols = [ind1,ind1,ind1,ind2,ind2,ind2];
            tarRows = [ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise,ind+2+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise,ind+2+s.SoundData.WhiteNoise];
        elseif ind > length(s.SoundData.UniqueFrequencies) - s.SoundData.WhiteNoise - 1
            tarCols = [ind1,ind1,ind1,ind2,ind2,ind2];
            tarRows = [ind-2+s.SoundData.WhiteNoise,ind-1+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise,ind-2+s.SoundData.WhiteNoise,ind-1+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise];
        end
        %now we need to pull the trial numbers from these targeting
        %systems!
        tarCount = 1;
        tarStore = [];
        for m = 1:length(tarCols)
            tarInds = find(trialMatrix(:,3) == uniqueDBs(tarCols(m)) & trialMatrix(:,2) == s.SoundData.UniqueFrequencies(tarRows(m)));
            tarStore(tarCount:tarCount + length(tarInds) - 1) = tarInds;
            tarCount = tarCount + length(tarInds);
        end
        tarStore = sort(tarStore);
        %now pull correct raster data!
        binSize = 0.0005;
        makeRast = functionBasicRaster(s.(desigName{j}).SpikeTimes,newSoundNon(tarStore),boundsLatHist);
        tarSpikes = makeRast(:,1);
%         findRast = ismember(s.(desigName{j}).AllRasters(:,2),tarStore);
%         tarSpikes = s.(desigName{j}).AllRasters(findRast,1);
        tarSpikes(tarSpikes < boundsLatHist(1)) = [];
        tarSpikes(tarSpikes > boundsLatHist(2)) = [];
        latFineHist(bigMasterInd + j - 1,:) = hist(tarSpikes,[boundsLatHist(1):binSize:boundsLatHist(2)])/binSize/length(tarStore);

        disp('About to do fine overall histogram')
        %generate a finer scale histogram across all tones
        makeRast = functionBasicRaster(s.(desigName{j}).SpikeTimes,newSoundNon(toneFinder),boundsLatHist);
        tarSpikes = makeRast(:,1);
%         findTars = find(ismember(s.(desigName{j}).AllRasters(:,2),toneFinder));
%         tarSpikes = s.(desigName{j}).AllRasters(findTars,1);
        tarSpikes(tarSpikes < boundsLatHist(1)) = [];
        tarSpikes(tarSpikes > boundsLatHist(2)) = [];
        fineHist(bigMasterInd + j - 1,:) = hist(tarSpikes,[boundsLatHist(1):binSize:boundsLatHist(2)])/binSize/length(toneFinder);
        
        %DO THE SAME FOR LASER
        tarCount = 1;
        tarStore = [];
        for m = 1:length(tarCols)
            tarInds = find(trialMatrixLaser(:,3) == uniqueDBs(tarCols(m)) & trialMatrixLaser(:,2) == s.SoundData.UniqueFrequencies(tarRows(m)));
            tarStore(tarCount:tarCount + length(tarInds) - 1) = tarInds;
            tarCount = tarCount + length(tarInds);
        end
        tarStore = sort(tarStore);
        %now pull correct raster data!
        binSize = 0.0005;
        makeRast = functionBasicRaster(s.(desigName{j}).SpikeTimes,newSoundLaser(tarStore),boundsLatHist);
        tarSpikes = makeRast(:,1);
%         findRast = ismember(s.(desigName{j}).AllRastersLaser(:,2),tarStore);
%         tarSpikes = s.(desigName{j}).AllRastersLaser(findRast,1);
        latFineHistLaser(bigMasterInd + j - 1,:) = hist(tarSpikes,[boundsLatHist(1):binSize:boundsLatHist(2)])/binSize/length(tarStore);

        disp('About to do fine overall histogram LASER')
        %generate a finer scale histogram across all tones
        makeRast = functionBasicRaster(s.(desigName{j}).SpikeTimes,newSoundLaser(toneFinderLaser),boundsLatHist);
        tarSpikes = makeRast(:,1);
%         findTars = find(ismember(s.(desigName{j}).AllRastersLaser(:,2),toneFinderLaser));
%         tarSpikes = s.(desigName{j}).AllRastersLaser(findTars,1);
        tarSpikes(tarSpikes < boundsLatHist(1)) = [];
        tarSpikes(tarSpikes > boundsLatHist(2)) = [];
%         tempFineRast = functionBasicRaster(s.(desigName{j}).SpikeTimes,trialMatrixLaser(toneFinderLaser,1),[-0.2 0.4]);
        
        fineHistLaser(bigMasterInd + j - 1,:) = hist(tarSpikes,[boundsLatHist(1):binSize:boundsLatHist(2)])/binSize/length(toneFinderLaser);
        
        nameStore{bigMasterInd + j -1} = targetFiles{i};
        unitStore{bigMasterInd + j -1} = desigName{j};
        tempBinStore(:,:,j) = s.(desigName{j}).BinDiff(1+s.SoundData.WhiteNoise:end,:,2);
        binValBigStore{bigMasterInd + j - 1} = s.(desigName{j}).BinDiff(1+s.SoundData.WhiteNoise:end,:,2);
        sigValBigStore{bigMasterInd + j - 1}  = s.(desigName{j}).BinSigVals(1+s.SoundData.WhiteNoise:end,:,2);
        binValBigStoreLaser{bigMasterInd + j - 1}  = s.(desigName{j}).BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2);
        sigValBigStoreLaser{bigMasterInd + j - 1}  = s.(desigName{j}).BinSigValsLaser(1+s.SoundData.WhiteNoise:end,:,2);
    end
    freqNumStore(bigMasterInd:bigMasterInd + numUnits - 1) = length(s.SoundData.UniqueFrequencies) - s.SoundData.WhiteNoise;
    bigDBStore(bigMasterInd:bigMasterInd + numUnits - 1) = dbStore;
    bigFreqStore(bigMasterInd:bigMasterInd + numUnits - 1) = freqStore;
    %now lets pull the correlation coefficients from tuning curves. 
    
    
    bigMasterInd = bigMasterInd + numUnits;
end

findMSNs = find(bigMaster(:,7) == 0);
findPVs = find(bigMaster(:,7) == 1);

% 

%% calculate widths with height based system
widthPercent = 0.5;

bigWidthHeightStore = [];
for i = 1:length(bigMaster)
    testValues = binValBigStore{i};
    dbStep(i) = 3/size(binValBigStore{i},1);
    for j = 1:size(testValues,2)
        [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues(:,j),widthPercent,3);
        bigWidthHeightStore(j,:,i) = widthVals;
        bigWidthMaxPosStore(j,i) = maxPos;
        bigWidthMaxValStore(j,i) = maxVal;
        bigWidthCutVal(j,i) = cutVal;
    end
end

%now lets just pull out the outward search, since this looks to be far more
%accurate
bigWidthSelWidth = bigWidthHeightStore(:,[1,3],:);

%repeat for laser trials
bigWidthHeightStoreLaser = [];
for i = 1:length(bigMaster)
    testValues = binValBigStoreLaser{i};
    for j = 1:size(testValues,2)
        [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues(:,j),widthPercent,3);
        bigWidthHeightStoreLaser(j,:,i) = widthVals;
        bigWidthMaxPosStoreLaser(j,i) = maxPos;
        bigWidthMaxValStoreLaser(j,i) = maxVal;
        bigWidthCutValLaser(j,i) = cutVal;
    end
end

%now lets just pull out the outward search, since this looks to be far more
%accurate
bigWidthSelWidthLaser = bigWidthHeightStoreLaser(:,[1,3],:);


%% Lets look at width values in normal trials

%now lets determine whether significant values exist at each DB level.
signStore= [];
signSig = [];
signSigStore = [];
sigCut = 0.05;
for i = 1:length(binValBigStore)
    signStore = sign(binValBigStore{i});
    signSig = signStore.* sigValBigStore{i};
    signSig(signSig <= 0) = NaN;
    signSigBina = min(signSig);
    signSigBina(min(signSig) < sigCut) = 10;
    signSigBina(signSigBina < 10) = 0;
    signSigBina(signSigBina ==10) = 1;
    findNan = isnan(signSigBina);
    signSigBina(findNan) = 0;
    signSigStore{i} = min(signSig);
    signSigBinaryStore{i} = signSigBina;
%     signSigBinaryStore(min(signSig) < sigCut,i) = 10;
%     signSigBinaryStore(signSigBinaryStore(:,i) < 10,i) = 0;
%     signSigBinaryStore(signSigBinaryStore(:,i) == 10,i) = 1;
end

%ok, now we have significant values, max points, and the detected width in
%each direction. Lets try and pull things out now?

sigWidthStore = cell(length(signSigBinaryStore),1);
sideWarnLow = cell(length(signSigBinaryStore),1);
sideWarnHi = cell(length(signSigBinaryStore),1);
for i = 1:length(binValBigStore)
    for j = 1:size(signSigBinaryStore{i},2)
        %first, check if there is a significant value
        if signSigBinaryStore{i}(j) == 1
            %now check for width values. If NaN on either side, determine
            %edge.
            %start with lower edge. 
            if isnan(bigWidthSelWidth(j,1,i))
                lowWidth = bigWidthMaxPosStore(j,i);
                sideWarnLow{i}(j) = 1;
            else
                lowWidth = bigWidthSelWidth(j,1,i);
            end
            if isnan(bigWidthSelWidth(j,2,i))
                hiWidth = numFreqUnit(i) - bigWidthMaxPosStore(j,i);
                sideWarnHi{i}(j) = 1;
            else
                hiWidth = bigWidthSelWidth(j,2,i);
            end
            sigWidthStore{i}(j) = (hiWidth + lowWidth)*dbStep(i);
        else
            sigWidthStore{i}(j) = 0;
        end
    end
end

%now lets just do some ugly stuff and look only at 70dB point. 
for i = 1:length(sigWidthStore)
    width70DB(i) = sigWidthStore{i}(end);
end


%% Do the same for laser trials

%now lets determine whether significant values exist at each DB level.
signStoreLaser= [];
signSig = [];
signSigStoreLaser = [];
sigCut = 0.05;
for i = 1:length(binValBigStoreLaser)
    signStoreLaser = sign(binValBigStoreLaser{i});
    signSig = signStoreLaser.* sigValBigStoreLaser{i};
    signSig(signSig <= 0) = NaN;
    signSigBina = min(signSig);
    signSigBina(min(signSig) < sigCut) = 10;
    signSigBina(signSigBina < 10) = 0;
    signSigBina(signSigBina ==10) = 1;
    findNan = isnan(signSigBina);
    signSigBina(findNan) = 0;
    signSigStoreLaser{i} = min(signSig);
    signSigBinaryStoreLaser{i} = signSigBina;
%     signSigBinaryStoreLaser(min(signSig) < sigCut,i) = 10;
%     signSigBinaryStoreLaser(signSigBinaryStoreLaser(:,i) < 10,i) = 0;
%     signSigBinaryStoreLaser(signSigBinaryStoreLaser(:,i) == 10,i) = 1;
end

%ok, now we have significant values, max points, and the detected width in
%each direction. Lets try and pull things out now?

sigWidthStoreLaser = cell(length(signSigBinaryStoreLaser),1);
sideWarnLowLaser = cell(length(signSigBinaryStoreLaser),1);
sideWarnHiLaser = cell(length(signSigBinaryStoreLaser),1);
for i = 1:length(binValBigStoreLaser)
    for j = 1:size(signSigBinaryStoreLaser{i},2)
        %first, check if there is a significant value
        if signSigBinaryStoreLaser{i}(j) == 1
            %now check for width values. If NaN on either side, determine
            %edge.
            %start with lower edge. 
            if isnan(bigWidthSelWidth(j,1,i))
                lowWidth = bigWidthMaxPosStoreLaser(j,i);
                sideWarnLowLaser{i}(j) = 1;
            else
                lowWidth = bigWidthSelWidth(j,1,i);
            end
            if isnan(bigWidthSelWidth(j,2,i))
                hiWidth = numFreqUnit(i) - bigWidthMaxPosStoreLaser(j,i);
                sideWarnHiLaser{i}(j) = 1;
            else
                hiWidth = bigWidthSelWidth(j,2,i);
            end
            sigWidthStoreLaser{i}(j) = (hiWidth + lowWidth)*dbStep(i);
        else
            sigWidthStoreLaser{i}(j) = 0;
        end
    end
end

%now lets just do some ugly stuff and look only at 70dB point. 
for i = 1:length(sigWidthStoreLaser)
    width70DBLaser(i) = sigWidthStoreLaser{i}(end);
end


%so just comparing all MSNs laser vs non, see a difference in width that is
%significant.
signrank(width70DB(findMSNs),width70DBLaser(findMSNs))
%this difference is a widening. 
mean(width70DB(findMSNs)) - mean(width70DBLaser(findMSNs))

%there is a larger mean effect on width in FSIs
mean(width70DB(findPVs)) - mean(width70DBLaser(findPVs))
%but no significance. 
signrank(width70DB(findPVs),width70DBLaser(findPVs))

%% find "good units"

%okay, thats significant, but what about selecting for good units?
sigCut = 0.05;
for i = 1:length(binValBigStoreLaser)
    %determine when sig and pos for baseline
    signStore = sign(binValBigStore{i});
    signSig = signStore.* sigValBigStore{i};
    numSig(i) = length(find(signSig < sigCut & signSig > 0));
    %determine when sig and pos for other
    signStoreLaser = sign(binValBigStoreLaser{i});
    signSigLaser = signStoreLaser.* sigValBigStoreLaser{i};
    numSigLaser(i) = length(find(signSigLaser < sigCut & signSigLaser > 0));
    %determine overlap matrix of laser and non-laser (for pulling stuff
    %like latency)
    signSig(signSig <= 0) = NaN;
    signSigLaser(signSigLaser <= 0) = NaN;
    tmpHold = signSig.*signSigLaser;
    tmpHold(~isnan(tmpHold)) = 1;
    overlapMatrix{i} = tmpHold;
end

%now lets go through and see which units had significant changes to size of
%response
for i = 1:length(binValBigStoreLaser)
    sigAmpStore(i,1) = signrank(reshape(binValBigStoreLaser{i},1,[]),reshape(binValBigStore{i},1,[]));
    testData = binValBigStoreLaser{i} - binValBigStore{i};
    changeDir(i) = sign(mean(mean(testData)));
    sigAmpStore(i,2) = signrank(reshape(testData,1,[]));
end



%also pull cells that have significant number of responses at baseline.
cutVal5 = 0.01;
cutValNum = 1;
for i = 1:length(binValBigStore)
    tmpData = sigValBigStore{i}; 
    tmpSign = sign(binValBigStore{i});
    tmpData = tmpData.*tmpSign;
    tmpData(tmpData <= 0) = 10;
    lengthSig(i) = length(find(tmpData <= cutVal5));
    altLengthSig(i) = length(find(sigValBigStore{i} <= 0.05));
end

for i = 1:length(binValBigStoreLaser)
    tmpData = sigValBigStoreLaser{i}; 
    tmpSign = sign(binValBigStoreLaser{i});
    tmpData = tmpData.*tmpSign;
    tmpData(tmpData <= 0) = 10;
    lengthSigLaser(i) = length(find(tmpData <= cutVal5));
end

%now lets determine the presence of significance at a particular dB level. 
cutVal5 = 0.01;
for i = 1:length(binValBigStoreLaser)
    if numDBUnit(i) == 3
        tmpSig = sigValBigStore{i};
    else
        tmpSig = sigValBigStore{i}(:,[1,5,9]);
    end
    tmpSig = min(tmpSig);
    tmpSig(tmpSig <= cutVal5) = 1;
    tmpSig(tmpSig < 1) = 0;
    tmpSig(isnan(tmpSig)) = 0;
    sigDBStore(:,i) = tmpSig;
end

%select units
tarCellsLight = find(sigAmpStore(:,1) < 0.05);
% tarCells = find(lengthSig > 5);
altCells = find(altLengthSig > 5);
tarCells = find(lengthSig > cutValNum | lengthSigLaser > cutValNum);
tarCells = intersect(tarCells,tarCellsLight);

tarPVs = intersect(tarCells,findPVs);
tarPVsDown = tarPVs;
tarPVsDown(changeDir(tarPVs)>0) = [];
tarMSNs = intersect(tarCells,findMSNs);
tarMSNsUp = tarMSNs;
tarMSNsUp(changeDir(tarMSNs) < 1) = [];
altPVs = intersect(altCells,findPVs);
altMSNs = intersect(altCells,findMSNs);

figure
axisLim = ceil(sqrt(length(tarMSNs)));
for i = 1:length(tarMSNs)
subplot(axisLim,axisLim,i)
plot(smooth(fineHist(tarMSNs(i),:),21))
end

% signrank(width70DB(tarMSNs),width70DBLaser(tarMSNs))

%% now try extracting from multiple DB bands. 

for i = 1:length(sigWidthStore)
    if numDBUnit(i) == 3
        widths(i,:) = sigWidthStore{i};
    else
        widths(i,:) = sigWidthStore{i}([1,5,9]);
    end
end

for i = 1:length(sigWidthStoreLaser)
    if numDBUnit(i) == 3
        widthsLaser(i,:) = sigWidthStoreLaser{i};
    else
        widthsLaser(i,:) = sigWidthStoreLaser{i}([1,5,9]);
    end
end

nanmean(widths(tarMSNs,1) - widthsLaser(tarMSNs,1))
signrank(widths(tarMSNs,1),widthsLaser(tarMSNs,1))
nanmean(widths(tarMSNs,2) - widthsLaser(tarMSNs,2))
signrank(widths(tarMSNs,2),widthsLaser(tarMSNs,2))
nanmean(widths(tarMSNs,3) - widthsLaser(tarMSNs,3))
signrank(widths(tarMSNs,3),widthsLaser(tarMSNs,3))

signrank(widths(tarPVs,1),widthsLaser(tarPVs,1))
signrank(widths(tarPVs,2),widthsLaser(tarPVs,2))
signrank(widths(tarPVs,3),widthsLaser(tarPVs,3))

%FSI effect significant. MSN effects significant if only selecting for
%response to light. 


%% Lets calculate the modulation index for the baseline period.

%want period from -0.3 to 0. This is 601:1201.

fineHistVect = [-0.6:0.0005:0.4];
modPV = -1*(mean(fineHist(findPVs,[600:1201])')- mean(fineHistLaser(findPVs,[600:1201])'))./(mean(fineHist(findPVs,[600:1201])')+ mean(fineHistLaser(findPVs,[600:1201])'));
modMSN = -1*(mean(fineHist(findMSNs,[600:1201])')- mean(fineHistLaser(findMSNs,[600:1201])'))./(mean(fineHist(findMSNs,[600:1201])')+ mean(fineHistLaser(findMSNs,[600:1201])'));
modTonePV = -1*(mean(fineHist(findPVs,[1201:1401])')- mean(fineHistLaser(findPVs,[1201:1401])'))./(mean(fineHist(findPVs,[1201:1401])')+ mean(fineHistLaser(findPVs,[1201:1401])'));
modToneMSN = -1*(mean(fineHist(findMSNs,[1201:1401])')- mean(fineHistLaser(findMSNs,[1201:1401])'))./(mean(fineHist(findMSNs,[1201:1401])')+ mean(fineHistLaser(findMSNs,[1201:1401])'));


hFig = figure;
subplot(2,1,1)
hist(modPV,[-1:0.1:1])
title('PV Modulation Index')
xlim([-1.05 1.05])
set(gca,'TickDir','out')
subplot(2,1,2)
hist(modMSN,[-1:0.1:1])
title('MSN Modulation Index')
xlim([-1.05 1.05])
set(gca,'TickDir','out')

spikeGraphName = 'modIndexBaseline';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
subplot(2,1,1)
hist(modTonePV,[-1:0.1:1])
title('PV Modulation Index')
xlim([-1.05 1.05])
set(gca,'TickDir','out')
subplot(2,1,2)
hist(modToneMSN,[-1:0.1:1])
title('MSN Modulation Index')
xlim([-1.05 1.05])
set(gca,'TickDir','out')

spikeGraphName = 'modIndexTone';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%plot combined. 

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
set(hFig, 'Position', [80 80 1000 400])
subplot(1,3,1)
hold on
for i = 1:length(findPVs)
    plot([1,2],[modPV(i),modTonePV(i)],'k')
end
errorbar([1,2],[mean(modPV) mean(modTonePV)],[std(modPV)/sqrt(length(modPV)), std(modTonePV)/sqrt(length(modPV))],'r')
xlim([0.5 2.5])
ylim([-1 1])
set(gca,'TickDir','out')

subplot(2,6,3)
hist(modPV,[-1:0.1:1])
title('PV Mod IndBase')
xlim([-1.05 1.05])
set(gca,'TickDir','out')
ylim([0 15])

subplot(2,6,9)
hist(modTonePV,[-1:0.1:1])
title('PV Mod IndTone')
xlim([-1.05 1.05])
set(gca,'TickDir','out')
ylim([0 15])

subplot(2,6,4)
hist(modMSN,[-1:0.1:1])
title('MSN Mod IndBase')
xlim([-1.05 1.05])
set(gca,'TickDir','out')
ylim([0 105])

subplot(2,6,10)
hist(modToneMSN,[-1:0.1:1])
title('MSN Mod IndTone')
xlim([-1.05 1.05])
set(gca,'TickDir','out')
ylim([0 105])

subplot(1,3,3)
hold on
for i = 1:length(findMSNs)
    plot([1,2],[modMSN(i),modToneMSN(i)],'k')
end
errorbar([1,2],[mean(modMSN) mean(modToneMSN)],[std(modMSN)/sqrt(length(modMSN)), std(modToneMSN)/sqrt(length(modMSN))],'r')
xlim([0.5 2.5])
ylim([-1 1])
set(gca,'TickDir','out')

spikeGraphName = 'OverallModIndexValues';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%plot scatter with best fit
[b,bintr,bintjm] = gmregress(modMSN,modToneMSN,0.05);
hFig = figure;
hold on
plot(modMSN,modToneMSN,'k.')
plot([-1 1],[-1 1],'r')
title('Scatter of Mod Index, X base Y tone')
set(gca,'TickDir','out')

spikeGraphName = 'modScatter';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

tester = modToneMSN- modMSN;


%% now lets plot scatters of width. 
for i = 1:3
    tarDB = i;
    finder = find(widths(tarMSNs,tarDB) == 0);
    find2 = find(widthsLaser(tarMSNs,tarDB) == 0);
    finder = unique([finder;find2]);
    tester = widths(tarMSNs,tarDB);
    test2 = widthsLaser(tarMSNs,tarDB);
    tester(finder) = [];
    test2(finder) = [];
    signrank(tester,test2)
    hFig = figure;
    plot(tester,test2,'r.')
    hold on
    plot([0 max(tester)],[0 max(tester)],'k')
    set(gca,'TickDir','out')
    spikeGraphName = strcat('ScatterWidthFor',num2str(i),'NoZero');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

for i = 1:3
    tarDB = i;
    
    tester = widths(tarMSNs,tarDB);
    test2 = widthsLaser(tarMSNs,tarDB);
    
    signrank(tester,test2)
    hFig = figure;
    plot(tester,test2,'r.')
    hold on
    plot([0 max(tester)],[0 max(tester)],'k')
    set(gca,'TickDir','out')
    spikeGraphName = strcat('ScatterWidthFor',num2str(i),'WithZero');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

%eliminating zero width units basically eliminates significance. 


%% Now lets try looking at aligning tuning curves for targeted units

%This uses a scaling factor of 20 to align all datapoints from multiple
%datasets.
scaleNum = 20;
widthPer = widthPercent;

for bigInd = 1:3
    tarDB = bigInd;

    bigStoreMSN = NaN(length(tarMSNs),41*scaleNum);
    bigStoreLaserMSN = NaN(length(tarMSNs),41*scaleNum);


    dbDecoder = [1,5,9];
    for i = 1:length(tarMSNs)
        %first, get the right DB stuff out
        if numDBUnit(tarMSNs(i)) == 3
            tmpDB = tarDB;
        else
            tmpDB = dbDecoder(tarDB);
        end
        %second, we need to figure out the proper increment.
        tmpInc = 3/(numFreqUnit(tarMSNs(i))-1);
        %pull curves
        tmp1 = binValBigStore{tarMSNs(i)}(:,tmpDB);
        tmp2 = binValBigStoreLaser{tarMSNs(i)}(:,tmpDB);
        %scale curves to match proper scaling.
        tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarMSNs(i))-1)];
        tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
        tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
        %find relative maxes. 
        [maxVal1 maxInd1] = max(tmp1Int);
        [maxVal2 maxInd2] = max(tmp2Int);
        [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
        bfSave(i,:) = [maxInd1,maxInd2];
        %divide by maximum value
        if maxVal1 > 0
            bigStoreMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;
            bigStoreLaserMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp2Int/maxVal1;
        end

    end

    hFig = figure;
    subplot(3,1,1)
    plot(nanmean(bigStoreMSN),'k')
    hold on
    plot(nanmean(bigStoreLaserMSN),'g')
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    set(gca,'TickDir','out')
    xlim([360 480])
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserMSN);
    test = sum(tester);
    plot(test,'g')
    xlim([360 480])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    subplot(3,1,3)
    hold on
    for i = 1:length(bfSave)
        plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
    end
    plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    spikeGraphName = strcat(['AverageTuningAlignedToBaseAtDB-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')

    % %try smoothing
    % tmpSmoothWind = 5;
    % figure
    % plot(smooth(nanmean(bigStoreMSN),tmpSmoothWind))
    % hold on
    % plot(smooth(nanmean(bigStoreLaserMSN),tmpSmoothWind),'r')
    % xlim([(21 - 5)*scaleNum,(21 + 5)*scaleNum])


    bigStoreMSN = NaN(length(tarMSNs),41*scaleNum);
    bigStoreLaserMSN = NaN(length(tarMSNs),41*scaleNum);

    % tarDB = 2;
    dbDecoder = [1,5,9];
    for i = 1:length(tarMSNs)
        %first, get the right DB stuff out
        if numDBUnit(tarMSNs(i)) == 3
            tmpDB = tarDB;
        else
            tmpDB = dbDecoder(tarDB);
        end
        %second, we need to figure out the proper increment.
        tmpInc = 3/(numFreqUnit(tarMSNs(i))-1);
        %pull curves
        tmp1 = binValBigStore{tarMSNs(i)}(:,tmpDB);
        tmp2 = binValBigStoreLaser{tarMSNs(i)}(:,tmpDB);
        %scale curves to match proper scaling.
        tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarMSNs(i))-1)];
        tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
        tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
        %find relative maxes. 
        [maxVal1 maxInd1] = max(tmp1Int);
        [maxVal2 maxInd2] = max(tmp2Int);
        [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
        bfSave(i,:) = [maxInd1,maxInd2];
        %divide by maximum value
        if maxVal1 > 0
            bigStoreMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;

        end
        if maxVal2 > 0
            bigStoreLaserMSN(i,21*scaleNum-maxInd2+1:21*scaleNum-maxInd2+length(tmp1Int)) = tmp2Int/maxVal2;
        end

    end
    
    hFig = figure;
    subplot(3,1,1)
    plot(nanmean(bigStoreMSN),'k')
    hold on
    plot(nanmean(bigStoreLaserMSN),'g')
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToEachAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    xlim([360 480])
    set(gca,'TickDir','out')
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserMSN);
    test = sum(tester);
    plot(test,'g')
    xlim([360 480])
    set(gca,'TickDir','out')
    title('Number of Data Points')
    subplot(3,1,3)
    hold on
    for i = 1:length(bfSave)
        plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
    end
    plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    
    spikeGraphName = strcat(['AverageTuningAlignedToEach-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end

%% Do the same for FSIs

for bigInd = 1:3
    tarDB = bigInd;

    bigStorePV = NaN(length(tarPVs),41*scaleNum);
    bigStoreLaserPV = NaN(length(tarPVs),41*scaleNum);


    dbDecoder = [1,5,9];
    for i = 1:length(tarPVs)
        %first, get the right DB stuff out
        if numDBUnit(tarPVs(i)) == 3
            tmpDB = tarDB;
        else
            tmpDB = dbDecoder(tarDB);
        end
        %second, we need to figure out the proper increment.
        tmpInc = 3/(numFreqUnit(tarPVs(i))-1);
        %pull curves
        tmp1 = binValBigStore{tarPVs(i)}(:,tmpDB);
        tmp2 = binValBigStoreLaser{tarPVs(i)}(:,tmpDB);
        %scale curves to match proper scaling.
        tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarPVs(i))-1)];
        tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
        tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
        %find relative maxes. 
        [maxVal1 maxInd1] = max(tmp1Int);
        [maxVal2 maxInd2] = max(tmp2Int);
        [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
        bfSave(i,:) = [maxInd1,maxInd2];
        %divide by maximum value
        if maxVal1 > 0
            bigStorePV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;
            bigStoreLaserPV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp2Int/maxVal1;
        end

    end

    hFig = figure;
    subplot(3,1,1)
    plot(nanmean(bigStorePV),'k')
    hold on
    plot(nanmean(bigStoreLaserPV),'g')
    tester = nanmean(bigStorePV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserPV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    set(gca,'TickDir','out')
    xlim([360 480])
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStorePV);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserPV);
    test = sum(tester);
    plot(test,'g')
    xlim([360 480])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    subplot(3,1,3)
    hold on
    for i = 1:length(bfSave)
        plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
    end
    plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    spikeGraphName = strcat(['AveragePVTuningAlignedToBaseAtDB-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')

    % %try smoothing
    % tmpSmoothWind = 5;
    % figure
    % plot(smooth(nanmean(bigStorePV),tmpSmoothWind))
    % hold on
    % plot(smooth(nanmean(bigStoreLaserPV),tmpSmoothWind),'r')
    % xlim([(21 - 5)*scaleNum,(21 + 5)*scaleNum])


    bigStorePV = NaN(length(tarPVs),41*scaleNum);
    bigStoreLaserPV = NaN(length(tarPVs),41*scaleNum);

    % tarDB = 2;
    dbDecoder = [1,5,9];
    for i = 1:length(tarPVs)
        %first, get the right DB stuff out
        if numDBUnit(tarPVs(i)) == 3
            tmpDB = tarDB;
        else
            tmpDB = dbDecoder(tarDB);
        end
        %second, we need to figure out the proper increment.
        tmpInc = 3/(numFreqUnit(tarPVs(i))-1);
        %pull curves
        tmp1 = binValBigStore{tarPVs(i)}(:,tmpDB);
        tmp2 = binValBigStoreLaser{tarPVs(i)}(:,tmpDB);
        %scale curves to match proper scaling.
        tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarPVs(i))-1)];
        tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
        tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
        %find relative maxes. 
        [maxVal1 maxInd1] = max(tmp1Int);
        [maxVal2 maxInd2] = max(tmp2Int);
        [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
        bfSave(i,:) = [maxInd1,maxInd2];
        %divide by maximum value
        if maxVal1 > 0
            bigStorePV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;

        end
        if maxVal2 > 0
            bigStoreLaserPV(i,21*scaleNum-maxInd2+1:21*scaleNum-maxInd2+length(tmp1Int)) = tmp2Int/maxVal2;
        end

    end
    
    hFig = figure;
    subplot(3,1,1)
    plot(nanmean(bigStorePV),'k')
    hold on
    plot(nanmean(bigStoreLaserPV),'g')
    tester = nanmean(bigStorePV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserPV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToEachAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    xlim([360 480])
    set(gca,'TickDir','out')
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStorePV);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserPV);
    test = sum(tester);
    plot(test,'g')
    xlim([360 480])
    set(gca,'TickDir','out')
    title('Number of Data Points')
    subplot(3,1,3)
    hold on
    for i = 1:length(bfSave)
        plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
    end
    plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    
    spikeGraphName = strcat(['AveragePVTuningAlignedToEach-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end

%% Now do MSNs but only include traces that have significant responses at the db level in question. 

bfDiffVect = [-30:2:30];
for bigInd = 1:3
    tarDB = bigInd;

    bigStoreMSN = NaN(length(tarMSNs),41*scaleNum);
    bigStoreLaserMSN = NaN(length(tarMSNs),41*scaleNum);
    bfSave = [];

    dbDecoder = [1,5,9];
    for i = 1:length(tarMSNs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(tarMSNs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(tarMSNs(i))-1);
            %pull curves
            tmp1 = binValBigStore{tarMSNs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{tarMSNs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarMSNs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStoreMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;
                bigStoreLaserMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp2Int/maxVal1;
            end
        end

    end

    hFig = figure;
    
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStoreMSN),'k','LineWidth',2)
    plot(nanmean(bigStoreMSN)+nanstd(bigStoreMSN)/sqrt(length(tarMSNs)),'k')
    plot(nanmean(bigStoreMSN)-nanstd(bigStoreMSN)/sqrt(length(tarMSNs)),'k')
    
    plot(nanmean(bigStoreLaserMSN),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserMSN)+nanstd(bigStoreLaserMSN)/sqrt(length(tarMSNs)),'g')
    plot(nanmean(bigStoreLaserMSN)-nanstd(bigStoreLaserMSN)/sqrt(length(tarMSNs)),'g')
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    set(gca,'TickDir','out')
    xlim([400 440])
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserMSN);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    spikeGraphName = strcat(['AverageMSNDBSelTuningAlignedToBaseAtDB-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')

    % %try smoothing
    % tmpSmoothWind = 5;
    % figure
    % plot(smooth(nanmean(bigStoreMSN),tmpSmoothWind))
    % hold on
    % plot(smooth(nanmean(bigStoreLaserMSN),tmpSmoothWind),'r')
    % xlim([(21 - 5)*scaleNum,(21 + 5)*scaleNum])


    bigStoreMSN = NaN(length(tarMSNs),41*scaleNum);
    bigStoreLaserMSN = NaN(length(tarMSNs),41*scaleNum);
    bfSave = [];
    % tarDB = 2;
    dbDecoder = [1,5,9];
    for i = 1:length(tarMSNs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(tarMSNs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(tarMSNs(i))-1);
            %pull curves
            tmp1 = binValBigStore{tarMSNs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{tarMSNs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarMSNs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStoreMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;

            end
            if maxVal2 > 0
                bigStoreLaserMSN(i,21*scaleNum-maxInd2+1:21*scaleNum-maxInd2+length(tmp1Int)) = tmp2Int/maxVal2;
            end
        end

    end
    
    hFig = figure;
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStoreMSN),'k','LineWidth',2)
    plot(nanmean(bigStoreMSN)+nanstd(bigStoreMSN)/sqrt(length(tarMSNs)),'k')
    plot(nanmean(bigStoreMSN)-nanstd(bigStoreMSN)/sqrt(length(tarMSNs)),'k')
    
    plot(nanmean(bigStoreLaserMSN),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserMSN)+nanstd(bigStoreLaserMSN)/sqrt(length(tarMSNs)),'g')
    plot(nanmean(bigStoreLaserMSN)-nanstd(bigStoreLaserMSN)/sqrt(length(tarMSNs)),'g')
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToEachAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    xlim([400 440])
    set(gca,'TickDir','out')
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserMSN);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    set(gca,'TickDir','out')
    title('Number of Data Points')
    
    subplot(3,1,3)
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     hold on
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    
    spikeGraphName = strcat(['AverageMSNDBSelTuningAlignedToEach-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end

%% Now apply this restricted analysis to FSIs. May not change much. 
for bigInd = 1:3
    tarDB = bigInd;

    bigStorePV = NaN(length(tarPVs),41*scaleNum);
    bigStoreLaserPV = NaN(length(tarPVs),41*scaleNum);
    bfSave = [];

    dbDecoder = [1,5,9];
    for i = 1:length(tarPVs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(tarPVs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(tarPVs(i))-1);
            %pull curves
            tmp1 = binValBigStore{tarPVs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{tarPVs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarPVs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStorePV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;
                bigStoreLaserPV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp2Int/maxVal1;
            end
        end

    end

    hFig = figure;
    
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStorePV),'k','LineWidth',2)
    plot(nanmean(bigStorePV)+nanstd(bigStorePV)/sqrt(length(tarPVs)),'k')
    plot(nanmean(bigStorePV)-nanstd(bigStorePV)/sqrt(length(tarPVs)),'k')
    
    plot(nanmean(bigStoreLaserPV),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserPV)+nanstd(bigStoreLaserPV)/sqrt(length(tarPVs)),'g')
    plot(nanmean(bigStoreLaserPV)-nanstd(bigStoreLaserPV)/sqrt(length(tarPVs)),'g')
    tester = nanmean(bigStorePV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserPV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    set(gca,'TickDir','out')
    xlim([400 440])
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStorePV);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserPV);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    spikeGraphName = strcat(['AveragePVDBSelTuningAlignedToBaseAtDB-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')

    % %try smoothing
    % tmpSmoothWind = 5;
    % figure
    % plot(smooth(nanmean(bigStorePV),tmpSmoothWind))
    % hold on
    % plot(smooth(nanmean(bigStoreLaserPV),tmpSmoothWind),'r')
    % xlim([(21 - 5)*scaleNum,(21 + 5)*scaleNum])


    bigStorePV = NaN(length(tarPVs),41*scaleNum);
    bigStoreLaserPV = NaN(length(tarPVs),41*scaleNum);
    bfSave = [];
    % tarDB = 2;
    dbDecoder = [1,5,9];
    for i = 1:length(tarPVs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(tarPVs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(tarPVs(i))-1);
            %pull curves
            tmp1 = binValBigStore{tarPVs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{tarPVs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(tarPVs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStorePV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;

            end
            if maxVal2 > 0
                bigStoreLaserPV(i,21*scaleNum-maxInd2+1:21*scaleNum-maxInd2+length(tmp1Int)) = tmp2Int/maxVal2;
            end
        end

    end
    
    hFig = figure;
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStorePV),'k','LineWidth',2)
    plot(nanmean(bigStorePV)+nanstd(bigStorePV)/sqrt(length(tarPVs)),'k')
    plot(nanmean(bigStorePV)-nanstd(bigStorePV)/sqrt(length(tarPVs)),'k')
    
    plot(nanmean(bigStoreLaserPV),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserPV)+nanstd(bigStoreLaserPV)/sqrt(length(tarPVs)),'g')
    plot(nanmean(bigStoreLaserPV)-nanstd(bigStoreLaserPV)/sqrt(length(tarPVs)),'g')
    tester = nanmean(bigStorePV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserPV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToEachAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    xlim([400 440])
    set(gca,'TickDir','out')
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStorePV);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserPV);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    set(gca,'TickDir','out')
    title('Number of Data Points')
    
    subplot(3,1,3)
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     hold on
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    
    spikeGraphName = strcat(['AveragePVDBSelTuningAlignedToEach-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end

%% Now lets repeat with alternative unit selection
%% Now do MSNs but only include traces that have significant responses at the db level in question. 

bfDiffVect = [-30:2:30];
for bigInd = 1:3
    tarDB = bigInd;

    bigStoreMSN = NaN(length(altMSNs),41*scaleNum);
    bigStoreLaserMSN = NaN(length(altMSNs),41*scaleNum);
    bfSave = [];

    dbDecoder = [1,5,9];
    for i = 1:length(altMSNs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(altMSNs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(altMSNs(i))-1);
            %pull curves
            tmp1 = binValBigStore{altMSNs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{altMSNs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(altMSNs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStoreMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;
                bigStoreLaserMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp2Int/maxVal1;
            end
        end

    end

    hFig = figure;
    
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStoreMSN),'k','LineWidth',2)
    plot(nanmean(bigStoreMSN)+nanstd(bigStoreMSN)/sqrt(length(altMSNs)),'k')
    plot(nanmean(bigStoreMSN)-nanstd(bigStoreMSN)/sqrt(length(altMSNs)),'k')
    
    plot(nanmean(bigStoreLaserMSN),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserMSN)+nanstd(bigStoreLaserMSN)/sqrt(length(altMSNs)),'g')
    plot(nanmean(bigStoreLaserMSN)-nanstd(bigStoreLaserMSN)/sqrt(length(altMSNs)),'g')
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    set(gca,'TickDir','out')
    xlim([400 440])
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserMSN);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    spikeGraphName = strcat(['AltMSNAverageMSNDBSelTuningAlignedToBaseAtDB-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')

    % %try smoothing
    % tmpSmoothWind = 5;
    % figure
    % plot(smooth(nanmean(bigStoreMSN),tmpSmoothWind))
    % hold on
    % plot(smooth(nanmean(bigStoreLaserMSN),tmpSmoothWind),'r')
    % xlim([(21 - 5)*scaleNum,(21 + 5)*scaleNum])


    bigStoreMSN = NaN(length(altMSNs),41*scaleNum);
    bigStoreLaserMSN = NaN(length(altMSNs),41*scaleNum);
    bfSave = [];
    % tarDB = 2;
    dbDecoder = [1,5,9];
    for i = 1:length(altMSNs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(altMSNs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(altMSNs(i))-1);
            %pull curves
            tmp1 = binValBigStore{altMSNs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{altMSNs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(altMSNs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStoreMSN(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;

            end
            if maxVal2 > 0
                bigStoreLaserMSN(i,21*scaleNum-maxInd2+1:21*scaleNum-maxInd2+length(tmp1Int)) = tmp2Int/maxVal2;
            end
        end

    end
    
    hFig = figure;
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStoreMSN),'k','LineWidth',2)
    plot(nanmean(bigStoreMSN)+nanstd(bigStoreMSN)/sqrt(length(altMSNs)),'k')
    plot(nanmean(bigStoreMSN)-nanstd(bigStoreMSN)/sqrt(length(altMSNs)),'k')
    
    plot(nanmean(bigStoreLaserMSN),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserMSN)+nanstd(bigStoreLaserMSN)/sqrt(length(altMSNs)),'g')
    plot(nanmean(bigStoreLaserMSN)-nanstd(bigStoreLaserMSN)/sqrt(length(altMSNs)),'g')
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserMSN);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToEachAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    xlim([400 440])
    set(gca,'TickDir','out')
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserMSN);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    set(gca,'TickDir','out')
    title('Number of Data Points')
    
    subplot(3,1,3)
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     hold on
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    
    spikeGraphName = strcat(['AltMSNAverageMSNDBSelTuningAlignedToEach-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end

%% Now apply this restricted analysis to FSIs. May not change much. 
for bigInd = 1:3
    tarDB = bigInd;

    bigStorePV = NaN(length(altPVs),41*scaleNum);
    bigStoreLaserPV = NaN(length(altPVs),41*scaleNum);
    bfSave = [];

    dbDecoder = [1,5,9];
    for i = 1:length(altPVs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(altPVs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(altPVs(i))-1);
            %pull curves
            tmp1 = binValBigStore{altPVs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{altPVs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(altPVs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStorePV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;
                bigStoreLaserPV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp2Int/maxVal1;
            end
        end

    end

    hFig = figure;
    
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStorePV),'k','LineWidth',2)
    plot(nanmean(bigStorePV)+nanstd(bigStorePV)/sqrt(length(altPVs)),'k')
    plot(nanmean(bigStorePV)-nanstd(bigStorePV)/sqrt(length(altPVs)),'k')
    
    plot(nanmean(bigStoreLaserPV),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserPV)+nanstd(bigStoreLaserPV)/sqrt(length(altPVs)),'g')
    plot(nanmean(bigStoreLaserPV)-nanstd(bigStoreLaserPV)/sqrt(length(altPVs)),'g')
    tester = nanmean(bigStorePV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserPV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    set(gca,'TickDir','out')
    xlim([400 440])
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStorePV);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserPV);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    spikeGraphName = strcat(['AltPVAveragePVDBSelTuningAlignedToBaseAtDB-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')

    % %try smoothing
    % tmpSmoothWind = 5;
    % figure
    % plot(smooth(nanmean(bigStorePV),tmpSmoothWind))
    % hold on
    % plot(smooth(nanmean(bigStoreLaserPV),tmpSmoothWind),'r')
    % xlim([(21 - 5)*scaleNum,(21 + 5)*scaleNum])


    bigStorePV = NaN(length(altPVs),41*scaleNum);
    bigStoreLaserPV = NaN(length(altPVs),41*scaleNum);
    bfSave = [];
    % tarDB = 2;
    dbDecoder = [1,5,9];
    for i = 1:length(altPVs)
        if sigDBStore(bigInd,i) == 1
            %first, get the right DB stuff out
            if numDBUnit(altPVs(i)) == 3
                tmpDB = tarDB;
            else
                tmpDB = dbDecoder(tarDB);
            end
            %second, we need to figure out the proper increment.
            tmpInc = 3/(numFreqUnit(altPVs(i))-1);
            %pull curves
            tmp1 = binValBigStore{altPVs(i)}(:,tmpDB);
            tmp2 = binValBigStoreLaser{altPVs(i)}(:,tmpDB);
            %scale curves to match proper scaling.
            tmpScale = [0:tmpInc*scaleNum:(tmpInc*scaleNum)*(numFreqUnit(altPVs(i))-1)];
            tmp1Int = interp1(tmpScale,tmp1,[1:(tmpScale(end))]);
            tmp2Int = interp1(tmpScale,tmp2,[1:(tmpScale(end))]);
            %find relative maxes. 
            [maxVal1 maxInd1] = max(tmp1Int);
            [maxVal2 maxInd2] = max(tmp2Int);
            [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
            bfSave(i,:) = [maxInd1,maxInd2];
            %divide by maximum value
            if maxVal1 > 0
                bigStorePV(i,21*scaleNum-maxInd1+1:21*scaleNum-maxInd1+length(tmp1Int)) = tmp1Int/maxVal1;

            end
            if maxVal2 > 0
                bigStoreLaserPV(i,21*scaleNum-maxInd2+1:21*scaleNum-maxInd2+length(tmp1Int)) = tmp2Int/maxVal2;
            end
        end

    end
    
    hFig = figure;
    subplot(3,1,1)
    hold on
    plot(nanmean(bigStorePV),'k','LineWidth',2)
    plot(nanmean(bigStorePV)+nanstd(bigStorePV)/sqrt(length(altPVs)),'k')
    plot(nanmean(bigStorePV)-nanstd(bigStorePV)/sqrt(length(altPVs)),'k')
    
    plot(nanmean(bigStoreLaserPV),'g','LineWidth',2)
    plot(nanmean(bigStoreLaserPV)+nanstd(bigStoreLaserPV)/sqrt(length(altPVs)),'g')
    plot(nanmean(bigStoreLaserPV)-nanstd(bigStoreLaserPV)/sqrt(length(altPVs)),'g')
    tester = nanmean(bigStorePV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthNorm = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreLaserPV);
    find1 = find(tester(1:420) <= widthPer,1,'last');
    widthLaser = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToEachAtDB-',num2str(tarDB),'WidthBase:',num2str(widthNorm),'WidthLaser:',num2str(widthLaser)])
    xlim([400 440])
    set(gca,'TickDir','out')
    
    subplot(3,1,2)
    hold on
    tester = ~isnan(bigStorePV);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreLaserPV);
    test = sum(tester);
    plot(test,'g')
    xlim([400 440])
    set(gca,'TickDir','out')
    title('Number of Data Points')
    
    subplot(3,1,3)
    
    subplot(3,1,3)
    hold on
    tmpZero = sum(bfSave');
    tmpZero = find(tmpZero == 0);
    tmpBF = bfSave;
    tmpBF(tmpZero,:) = [];
    tmpBFDiff = tmpBF(:,1) - tmpBF(:,2);
    hist(tmpBFDiff,bfDiffVect)
%     hold on
%     for i = 1:length(bfSave)
%         plot([1,2],[bfSave(i,1),bfSave(i,2)],'k')
%     end
%     plot([1,2],mean(bfSave),'r','LineWidth',2)
    title(['signrank-',num2str(signrank(bfSave(:,1),bfSave(:,2))),'~meandiff-',num2str(mean(bfSave(:,1) - mean(bfSave(:,2)))),'-std-',num2str(std(bfSave(:,1) - bfSave(:,2)))])
    
    spikeGraphName = strcat(['AltPVAveragePVDBSelTuningAlignedToEach-',num2str(tarDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end


%% Now lets plot out example units

load('180918_ML180618E_R_AudStr_3300_3mWPVHaloTuningWhiteAltLaserFullTuningAnalysis')

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
set(hFig, 'Position', [80 80 1000 400])

subplot(4,8,17)
%find max waveform
waves = s.nt10cluster3.AverageWaveForms;
[maxVal maxInd] = max(max(waves));
plot(waves(:,maxInd))
set(gca,'TickDir','out')

%plot overall histograms
subplot(2,4,1)
hold on
plot(s.nt10cluster3.HistBinVector,s.nt10cluster3.AllHistograms,'k')
plot(s.nt10cluster3.HistBinVector,s.nt10cluster3.AllHistogramsLaser,'g')
set(gca,'TickDir','out')

%plot scatter of binned values
subplot(2,4,2)
hold on
plot(s.nt10cluster3.BinDiff(1+s.SoundData.WhiteNoise:end,:,2),s.nt10cluster3.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2),'k.')
absMin = min(min(min(s.nt10cluster3.BinDiff(1+s.SoundData.WhiteNoise:end,:,2)),min(s.nt10cluster3.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2))));
absMax = max(max(max(s.nt10cluster3.BinDiff(1+s.SoundData.WhiteNoise:end,:,2)),max(s.nt10cluster3.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2))));
plot([absMin absMax],[absMin absMax],'r')
xlim([absMin absMax])
ylim([absMin absMax])
set(gca,'TickDir','out')

%plot wall of hist
subplot(1,2,2)
hold on

wallHolder = [];
wallSpacer = 10;
wallSmoothWind = 7;
% counter = 1;
for i = 1:3
    counter = 1;
    for j = 1:16
        wallHolder(counter:counter + 40,i) = smooth(s.nt10cluster3.FreqDBHistograms(j+s.SoundData.WhiteNoise,i,81:121),wallSmoothWind);
        counter = counter + 40;
        wallHolder(counter:counter + length(wallSpacer)-1,i) = NaN;
        counter = counter + length(wallSpacer);
    end
end

wallHolderLaser = [];
wallSpacer = 10;
wallSmoothWind = 7;
% counter = 1;
for i = 1:3
    counter = 1;
    for j = 1:16
        wallHolderLaser(counter:counter + 40,i) = smooth(s.nt10cluster3.FreqDBHistogramsLaser(j+s.SoundData.WhiteNoise,i,81:121),wallSmoothWind);
        counter = counter + 40;
        wallHolderLaser(counter:counter + length(wallSpacer)-1,i) = NaN;
        counter = counter + length(wallSpacer);
    end
end

wallMax = max(max(wallHolder));
wallMaxLaser = max(max(wallHolderLaser));

for i = 1:3
    plot(wallHolder(:,i)/wallMax + i - 1,'k')
    plot(wallHolderLaser(:,i)/wallMax + i - 1,'g')
end
set(gca,'TickDir','out')

title(num2str(wallMax))



spikeGraphName = 'FSIExampleWallHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')





hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
set(hFig, 'Position', [80 80 1000 400])

subplot(4,8,17)
%find max waveform
waves = s.nt13cluster2.AverageWaveForms;
[maxVal maxInd] = max(max(waves));
plot(waves(:,maxInd))
set(gca,'TickDir','out')
%plot overall histograms
subplot(2,4,1)
hold on
plot(s.nt13cluster2.HistBinVector,s.nt13cluster2.AllHistograms,'k')
plot(s.nt13cluster2.HistBinVector,s.nt13cluster2.AllHistogramsLaser,'g')
set(gca,'TickDir','out')
%plot scatter of binned values
subplot(2,4,2)
hold on
plot(s.nt13cluster2.BinDiff(1+s.SoundData.WhiteNoise:end,:,2),s.nt13cluster2.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2),'k.')
absMin = min(min(min(s.nt13cluster2.BinDiff(1+s.SoundData.WhiteNoise:end,:,2)),min(s.nt13cluster2.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2))));
absMax = max(max(max(s.nt13cluster2.BinDiff(1+s.SoundData.WhiteNoise:end,:,2)),max(s.nt13cluster2.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2))));
plot([absMin absMax],[absMin absMax],'r')
xlim([absMin absMax])
ylim([absMin absMax])
set(gca,'TickDir','out')
%plot wall of hist
subplot(1,2,2)
hold on

wallHolder = [];
wallSpacer = 10;
wallSmoothWind = 7;
% counter = 1;
for i = 1:3
    counter = 1;
    for j = 1:16
        wallHolder(counter:counter + 40,i) = smooth(s.nt13cluster2.FreqDBHistograms(j+s.SoundData.WhiteNoise,i,81:121),wallSmoothWind);
        counter = counter + 40;
        wallHolder(counter:counter + length(wallSpacer)-1,i) = NaN;
        counter = counter + length(wallSpacer);
    end
end

wallHolderLaser = [];
wallSpacer = 10;
wallSmoothWind = 7;
% counter = 1;
for i = 1:3
    counter = 1;
    for j = 1:16
        wallHolderLaser(counter:counter + 40,i) = smooth(s.nt13cluster2.FreqDBHistogramsLaser(j+s.SoundData.WhiteNoise,i,81:121),wallSmoothWind);
        counter = counter + 40;
        wallHolderLaser(counter:counter + length(wallSpacer)-1,i) = NaN;
        counter = counter + length(wallSpacer);
    end
end

wallMax = max(max(wallHolder));
wallMaxLaser = max(max(wallHolderLaser));

for i = 1:3
    plot(wallHolder(:,i)/wallMaxLaser + i - 1,'k')
    plot(wallHolderLaser(:,i)/wallMaxLaser + i - 1,'g')
end
set(gca,'TickDir','out')
title(num2str(wallMaxLaser))



spikeGraphName = 'MSNAdditiveExampleWallHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')







hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
set(hFig, 'Position', [80 80 1000 400])

subplot(4,8,17)
%find max waveform
waves = s.nt13cluster1.AverageWaveForms;
[maxVal maxInd] = max(max(waves));
plot(waves(:,maxInd))
set(gca,'TickDir','out')
%plot overall histograms
subplot(2,4,1)
hold on
plot(s.nt13cluster1.HistBinVector,s.nt13cluster1.AllHistograms,'k')
plot(s.nt13cluster1.HistBinVector,s.nt13cluster1.AllHistogramsLaser,'g')
set(gca,'TickDir','out')
%plot scatter of binned values
subplot(2,4,2)
hold on
plot(s.nt13cluster1.BinDiff(1+s.SoundData.WhiteNoise:end,:,2),s.nt13cluster1.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2),'k.')
absMin = min(min(min(s.nt13cluster1.BinDiff(1+s.SoundData.WhiteNoise:end,:,2)),min(s.nt13cluster1.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2))));
absMax = max(max(max(s.nt13cluster1.BinDiff(1+s.SoundData.WhiteNoise:end,:,2)),max(s.nt13cluster1.BinDiffLaser(1+s.SoundData.WhiteNoise:end,:,2))));
plot([absMin absMax],[absMin absMax],'r')
xlim([absMin absMax])
ylim([absMin absMax])
set(gca,'TickDir','out')
%plot wall of hist
subplot(1,2,2)
hold on

wallHolder = [];
wallSpacer = 10;
wallSmoothWind = 7;
% counter = 1;
for i = 1:3
    counter = 1;
    for j = 1:16
        wallHolder(counter:counter + 40,i) = smooth(s.nt13cluster1.FreqDBHistograms(j+s.SoundData.WhiteNoise,i,81:121),wallSmoothWind);
        counter = counter + 40;
        wallHolder(counter:counter + length(wallSpacer)-1,i) = NaN;
        counter = counter + length(wallSpacer);
    end
end

wallHolderLaser = [];
wallSpacer = 10;
wallSmoothWind = 7;
% counter = 1;
for i = 1:3
    counter = 1;
    for j = 1:16
        wallHolderLaser(counter:counter + 40,i) = smooth(s.nt13cluster1.FreqDBHistogramsLaser(j+s.SoundData.WhiteNoise,i,81:121),wallSmoothWind);
        counter = counter + 40;
        wallHolderLaser(counter:counter + length(wallSpacer)-1,i) = NaN;
        counter = counter + length(wallSpacer);
    end
end

wallMax = max(max(wallHolder));
wallMaxLaser = max(max(wallHolderLaser));

for i = 1:3
    plot(wallHolder(:,i)/wallMaxLaser + i - 1,'k')
    plot(wallHolderLaser(:,i)/wallMaxLaser + i - 1,'g')
end
set(gca,'TickDir','out')
title(num2str(wallMaxLaser))



spikeGraphName = 'MSNMult1ExampleWallHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%%Create bar plots of modulation levels

load('linMatCombined.mat')

maxMSN = sum(sum(linMatMSN));
maxFSI = sum(sum(linMatFSI));

hFig = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.01 0.01]);
set(hFig, 'Position', [80 80 700 700])

for i = 1:3
    for j = 1:3
        subplot(3,3,3*(i-1) + j)
        bar([1 2],[linMatFSI(i,j)/maxFSI,linMatMSN(i,j)/maxMSN])
        xlim([0.5 2.5])
        ylim([0 0.5])
        set(gca,'YTick',[])
        set(gca,'XTick',[])
    end
end


spikeGraphName = 'LinearSingleUnitPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')






