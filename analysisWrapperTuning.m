%This is code to do wrapper functions on analysisTuningWithWhite output
%data. Takes in s and masterData in, and is meant to output overall
%information about the population. This will be for data from baseline
%tuning recordings.


%% identify files, pull names, extraction via loop

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

%actually extract files.
for i = 1:numFiles
    %load the target file
    disp(strcat('Analyzing File:',num2str(i)))
    load(targetFiles{i})
    numUnits = size(masterData,1);
    numDBs = s.SoundData.NumDBs;
    pkTroughRatio = [];
    interpWaves = [];
    tempBinStore = [];
    tempSigStore = [];
    tempLatStore = [];
    tempMaxStore = [];
    tempWidthStore = [];
    tempHist = [];
    dbStore = [];
    freqStore = [];
    halfWidthTime = [];
    pkTroughTime = [];
    %generate a finder for non-white noise tones
    toneFinder = find(s.TrialMatrix(:,2) ~= 0);
%     toneFinder(toneFinder > length(s.SoundData.TrialMatrix)) = [];
    %now lets also generate a trial matrix for all tone times.
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
    %since non-white noise sessions will have funky tone DB values, round
    %up to the nearest ten.
    trialMatrix(:,3) = round(trialMatrix(:,3)/10)*10;
    uniqueDBs = unique(trialMatrix(:,3));
    %pull from masterData, store in overall.
    bigMaster(bigMasterInd:bigMasterInd + numUnits - 1,:) = masterData(:,1:12);
    %now pull overall designations of pos/neg/mix/no response
    [indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
    [indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');

    holder = masterData(:,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    respVect(:,bigMasterInd:bigMasterInd + numUnits - 1) = (holder(:,1) + holder(:,2))'; %note that here, -2 = neg, -1 = mix, 0 = no, 1 = pos
    %find the peak to trough value ratio in waveforms.
    numCells = length(s.DesignationName);
    desigName = s.DesignationName;
    for j = 1:numCells
        disp(desigName{j})
        %store RPVs
        rpvStore(bigMasterInd + j -1) = s.(desigName{j}).RPVPercent;
        rpvNumStore(bigMasterInd + j -1) = s.(desigName{j}).RPVNumber;
        %pull up cell average waveforms
        cellWaves = s.(desigName{j}).AverageWaveForms;
        %now pull up which one is biggest
        waveMax = max(cellWaves);
        [maxVal maxInd] = max(waveMax);
        %now generate finely sampled version
        chosenWave = cellWaves(:,maxInd);
        interpVect = [1:0.1:40];
        interpWave = interp1(1:40,chosenWave,interpVect,'spline');
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
        pkTroughTime(j) = troughInd;
        troughInd = troughInd + pkInd - 1;
        pkTroughRatio(j) = pkVal/troughVal;
        interpWaves(j,:) = interpWave;
        %now get half-width
        halfFirst = find(interpWave(1:pkInd) > pkVal/2,1,'first');
        halfSecond = find(interpWave(halfFirst:end) < pkVal/2,1,'first');
        halfWidthTime(j) = halfSecond;

        %pull out binned values for entire tone period, as well as
        %significance
        tempBinStore(:,:,j) = s.(desigName{j}).BinDiff(1+s.SoundData.WhiteNoise:end,end-5+1:end,2);
        tempSigStore(:,:,j) = s.(desigName{j}).BinSigVals(1+s.SoundData.WhiteNoise:end,end-5+1:end,2);
        tempLatStore(:,:,j) = s.(desigName{j}).LatencyMap(1+s.SoundData.WhiteNoise:end,end-5+1:end);
        try
            tempWidthStore(:,:,j) = s.WidthLatData(1+s.SoundData.WhiteNoise:end,end-5+1:end,j);
        catch
            tempWidthStore(:,:,j) = s.NonLaserOverall.WidthLatData(1+s.SoundData.WhiteNoise:end,end-5+1:end,j);
            
        end
        
        [row,col] = find(s.(desigName{j}).BinTone(1+s.SoundData.WhiteNoise:end,end-5+1:end) == max(max(s.(desigName{j}).BinTone(1+s.SoundData.WhiteNoise:end,end-5+1:end))));
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
        if ind >=2 & ind <= 15 %these are non-edge cases
            tarCols = [ind1,ind1,ind1,ind2,ind2,ind2];
            tarRows = [ind-1+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise,ind-1+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise];
        elseif ind < 2
            tarCols = [ind1,ind1,ind1,ind2,ind2,ind2];
            tarRows = [ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise,ind+2+s.SoundData.WhiteNoise,ind+s.SoundData.WhiteNoise,ind+1+s.SoundData.WhiteNoise,ind+2+s.SoundData.WhiteNoise];
        elseif ind > 15
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
        findRast = ismember(s.(desigName{j}).AllRasters(:,2),tarStore);
        tarSpikes = s.(desigName{j}).AllRasters(findRast,1);
        latFineHist(bigMasterInd + j - 1,:) = hist(tarSpikes,[-0.2:binSize:0.4])/binSize/length(tarStore);

        disp('About to do fine overall histogram')
        %generate a finer scale histogram across all tones
        tempFineRast = functionBasicRaster(s.(desigName{j}).SpikeTimes,s.TrialMatrix(toneFinder,1),[-0.2 0.4]);
        
        fineHist(bigMasterInd + j - 1,:) = hist(tempFineRast(:,1),[-0.2:binSize:0.4])/binSize/length(toneFinder);
        nameStore{bigMasterInd + j -1} = targetFiles{i};
        recNumStore(bigMasterInd+j-1) = i;
        unitStore{bigMasterInd + j -1} = desigName{j};
    end
    halfWidthTimeStore(bigMasterInd:bigMasterInd + numUnits - 1) = halfWidthTime;
    pkTroughTimeStore(bigMasterInd:bigMasterInd + numUnits - 1) = pkTroughTime;
    pkTroughRatioStore(bigMasterInd:bigMasterInd + numUnits - 1) = pkTroughRatio;
    interpWaveStore(bigMasterInd:bigMasterInd + numUnits - 1,:) = interpWaves;
    binValBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempBinStore;
    sigValBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempSigStore;
    latMapBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempLatStore;
    if size(tempWidthStore,1) < size(unique(s.TrialMatrix(:,2)),1) - s.SoundData.WhiteNoise
        sizeDiff = size(s.NonLaserOverall.WidthData,1) - size(tempWidthStore,1);
        tempWidthStore(end+sizeDiff,:,:) = zeros(sizeDiff,size(tempWidthStore,2),size(tempWidthStore,3));
    end
    widthLatStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempWidthStore;
%         bigMaxStore(bigMasterInd:bigMasterInd + numUnits - 1) =
%         tempMaxStore;
%     bigHistStore(:,bigMasterInd:bigMasterInd + numUnits - 1) = tempHist;
    bigDBStore(bigMasterInd:bigMasterInd + numUnits - 1) = dbStore;
    bigFreqStore(bigMasterInd:bigMasterInd + numUnits - 1) = freqStore;
    recStore(bigMasterInd:bigMasterInd + numUnits - 1) = i;

    %instead, lets just pull from posWidths and negWidths
    try
        intFastPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.PosWidths((end-5+1:end),:,1));
        intFastNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NegWidths((end-5+1:end),:,1));
        intSlowPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.PosWidths((end-5+1:end),:,3));
        intSlowNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NegWidths((end-5+1:end),:,3));
        %store positive tuning widths
        widthStore(:,bigMasterInd:bigMasterInd + numUnits - 1,:) = s.PosWidths(end-5+1:end,:,:);
    catch
        intFastPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NonLaserOverall.PosWidths((end-5+1:end),:,1));
        intFastNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NonLaserOverall.NegWidths((end-5+1:end),:,1));
        intSlowPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NonLaserOverall.PosWidths((end-5+1:end),:,3));
        intSlowNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NonLaserOverall.NegWidths((end-5+1:end),:,3));
        %store positive tuning widths
        widthStore(:,bigMasterInd:bigMasterInd + numUnits - 1,:) = s.NonLaserOverall.PosWidths(end-5+1:end,:,:);
    end
    
    %store BFs
    bfStore(bigMasterInd:bigMasterInd + numUnits - 1) = masterData(:,12);
    
    %now lets try to get corrcoef from binDiff values (tempBinStore)
    shuffBinDiffCorrVals = [];
    binDiffCorrVals = [];
    binDiffCorrValSig = [];
    for j = 1:numUnits
        val1 = reshape(squeeze(tempBinStore(:,:,j)),1,[]);
        for k = 1:numUnits
            val2 = reshape(squeeze(tempBinStore(:,:,k)),1,[]);
            tempCorrVal = corrcoef(val1,val2);
            binDiffCorrVals(j,k) = tempCorrVal(2);
            %now do shuffled versions
            numShuff = 1000;
            for l = 1:numShuff
                shuffVal = randperm(length(val2));
                shuffVal = val2(shuffVal);
                tempCorrVal = corrcoef(val1,shuffVal);
                shuffBinDiffCorrVals(l) = tempCorrVal(2);
            end
            %calculate percentile
            prctileVals = prctile(shuffBinDiffCorrVals,[0.5 99.5]);
            if binDiffCorrVals(j,k) < prctileVals(1) || binDiffCorrVals(j,k) > prctileVals(2)
                binDiffCorrValSig(j,k) = 1;
            else
                binDiffCorrValSig(j,k) = 0;
            end
            shuffBinDiffCorrVals = [];
            %also calculate distances between units
            distVal(j,k) = sqrt(((floor(masterData(j,1)) - floor(masterData(k,1)))*25)^2 + ((masterData(j,2) - masterData(k,2))*250)^2);
        end
    end
    findPVs = find(masterData(:,5) < 0.0004 & masterData(:,6) > 1.1);
    findMSNs = find(masterData(:,5) > 0.0005 & masterData(:,6) > 1.1); 
    findCHATs = find(masterData(:,6) < 1.1);
    corrCoefData.CorrCoefs = binDiffCorrVals;
    corrCoefData.CorrCoefSig = binDiffCorrValSig;
    corrCoefData.CellType = NaN(numUnits,1);
    corrCoefData.CellType(findPVs) = 1;
    corrCoefData.CellType(findMSNs) = 0;
    corrCoefData.CellType(findCHATs) = 4;
    corrCoefData.Distance = distVal;
    bigCorrStore{i} = corrCoefData;
    corrCoefData = [];
    distVal = [];
    
    
    bigMasterInd = bigMasterInd + numUnits;
        
end

%% determine unit balance for each recording

for i = 1:numFiles
    %find targets
    tarFind = find(recNumStore == i);
    %now find number of MSNs and FSIs
    cellRecStore(i,1) = length(find(bigMaster(tarFind,7) == 0));
    cellRecStore(i,2) = length(find(bigMaster(tarFind,7) == 1));
end

%% now lets march through the correlation coefficient data and see what we pull out. First things first, I want to remove all the excess values, since I dont want double counting. 
%lets define interactions! There are three cell types, and therefore six
%interactions: FSI-MSN, FSI-ChAT, FSI-FSI, MSN-ChAT, MSN-MSN, ChAT-ChAT. 

%lets number these. I want the most common first, going up into the less
%common. so:

% MSN-MSN: 0
% MSN-FSI: 1
% FSI-FSI: 2
% MSN-ChAT: 4
% FSI-ChAT: 5
% ChAT-ChAT: 8

%First, lets use this as a platform to determine the number of possible FSI
%-MSN pairs I can sample within a reasonable distance of each other
tarCorrs = [2,9,10,13,15,16];

corrCount = 1;
cellCount = 1;
for i = 1:length(tarCorrs)
    dataset = bigCorrStore{tarCorrs(i)};
    %now lets try and linearize this dataset. 
    numUnits = length(dataset.CellType);
    for j = 1:numUnits
        corrValStoreRestrict(corrCount:corrCount + numUnits - j - 1) = dataset.CorrCoefs(j,j+1:end);
        distValStoreRestrict(corrCount:corrCount + numUnits - j - 1) = dataset.Distance(j,j+1:end);
        interTypeStoreRestrict(corrCount:corrCount + numUnits - j - 1) = dataset.CellType(j)+ dataset.CellType(j+1:end);
        corrValSigStoreRestrict(corrCount:corrCount + numUnits - j - 1) = dataset.CorrCoefSig(j,j+1:end);
        cellTrackerStoreRestrict(corrCount:corrCount + numUnits - j - 1,1) = i;
        cellTrackerStoreRestrict(corrCount:corrCount + numUnits - j - 1,2) = j+cellCount-1;
        cellTrackerStoreRestrict(corrCount:corrCount + numUnits - j - 1,3) = cellCount + (j+1:numUnits) - 1;
        cellTrackerStoreRestrict(corrCount:corrCount + numUnits - j - 1,4) = dataset.CellType(j);
        cellTrackerStoreRestrict(corrCount:corrCount + numUnits - j - 1,5) = dataset.CellType(j+1:numUnits);
        cellTrackerStoreRestrict(corrCount:corrCount + numUnits - j - 1,6) = bigMaster(j+cellCount-1,2);
        cellTrackerStoreRestrict(corrCount:corrCount + numUnits - j - 1,7) = bigMaster(cellCount + (j+1:numUnits) - 1,2);
        
        corrCount = corrCount + numUnits - j;
    end
    cellCount = cellCount + numUnits;
end

%now lets clean up the dataset. Remove all NaNs. 
corrValStoreRestrict(isnan(interTypeStoreRestrict)) = [];
distValStoreRestrict(isnan(interTypeStoreRestrict)) = [];
cellTrackerStoreRestrict(isnan(interTypeStoreRestrict),:) = [];
corrValSigStoreRestrict(isnan(interTypeStoreRestrict)) = [];
interTypeStoreRestrict(isnan(interTypeStoreRestrict)) = [];

%MSN vs FSI, close
selDist = 10000;
selInterType = 1;


findTarCorr = intersect(find(distValStoreRestrict <= selDist),find(interTypeStoreRestrict == selInterType));


corrCount = 1;
cellCount = 1;
for i = 1:length(bigCorrStore)
    dataset = bigCorrStore{i};
    %now lets try and linearize this dataset. 
    numUnits = length(dataset.CellType);
    for j = 1:numUnits
        corrValStore(corrCount:corrCount + numUnits - j - 1) = dataset.CorrCoefs(j,j+1:end);
        distValStore(corrCount:corrCount + numUnits - j - 1) = dataset.Distance(j,j+1:end);
        interTypeStore(corrCount:corrCount + numUnits - j - 1) = dataset.CellType(j)+ dataset.CellType(j+1:end);
        corrValSigStore(corrCount:corrCount + numUnits - j - 1) = dataset.CorrCoefSig(j,j+1:end);
        cellTrackerStore(corrCount:corrCount + numUnits - j - 1,1) = i;
        cellTrackerStore(corrCount:corrCount + numUnits - j - 1,2) = j+cellCount-1;
        cellTrackerStore(corrCount:corrCount + numUnits - j - 1,3) = cellCount + (j+1:numUnits) - 1;
        cellTrackerStore(corrCount:corrCount + numUnits - j - 1,4) = dataset.CellType(j);
        cellTrackerStore(corrCount:corrCount + numUnits - j - 1,5) = dataset.CellType(j+1:numUnits);
        cellTrackerStore(corrCount:corrCount + numUnits - j - 1,6) = bigMaster(j+cellCount-1,2);
        cellTrackerStore(corrCount:corrCount + numUnits - j - 1,7) = bigMaster(cellCount + (j+1:numUnits) - 1,2);
        
        corrCount = corrCount + numUnits - j;
    end
    cellCount = cellCount + numUnits;
end

%now lets clean up the dataset. Remove all NaNs. 
corrValStore(isnan(interTypeStore)) = [];
distValStore(isnan(interTypeStore)) = [];
cellTrackerStore(isnan(interTypeStore),:) = [];
corrValSigStore(isnan(interTypeStore)) = [];
interTypeStore(isnan(interTypeStore)) = [];

%lets look in general at within neuron corrcoefs. 
corrHistVect = [-1:0.05:1];

%just MSNs, close
selDist = 150;
selInterType = 0;

findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
hFig = figure;
corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValmsnmsnSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValmsnmsn,'w')
bar(corrHistVect,corrValmsnmsnSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])
set(gca,'TickDir','out');
% corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);
xlabel('Correlation Coefficient')
ylabel('Number of Pairwise Comparisons')


spikeGraphName = 'TuningCorrelationMSN';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%just FSIs, close
selDist = 150;
selInterType = 2;

findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
hFig = figure;
corrValfsifsi = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValfsifsiSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValfsifsi,'w')
bar(corrHistVect,corrValfsifsiSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])
set(gca,'TickDir','out');
xlabel('Correlation Coefficient')
ylabel('Number of Pairwise Comparisons')


spikeGraphName = 'TuningCorrelationFSI';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%MSN vs FSI, close
selDist = 150;
selInterType = 1;


findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
hFig = figure;
corrValfsimsn = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValfsimsnSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValfsimsn,'w')
bar(corrHistVect,corrValfsimsnSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])
set(gca,'TickDir','out');
xlabel('Correlation Coefficient')
ylabel('Number of Pairwise Comparisons')


spikeGraphName = 'TuningCorrelationFSIMSN';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%NOW DO THESE ALL FAR
%just MSNs, close
selDist = 1200;
selInterType = 0;

findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
hFig = figure;
corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValmsnmsnSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValmsnmsn,'w')
bar(corrHistVect,corrValmsnmsnSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])
set(gca,'TickDir','out');
% corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);
xlabel('Correlation Coefficient')
ylabel('Number of Pairwise Comparisons')


spikeGraphName = 'allDistTuningCorrelationMSN';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%just FSIs, close
selDist = 1200;
selInterType = 2;

findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
hFig = figure;
corrValfsifsi = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValfsifsiSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValfsifsi,'w')
bar(corrHistVect,corrValfsifsiSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])
set(gca,'TickDir','out');
xlabel('Correlation Coefficient')
ylabel('Number of Pairwise Comparisons')


spikeGraphName = 'allDistTuningCorrelationFSI';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%MSN vs FSI, close
selDist = 1200;
selInterType = 1;


findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
hFig = figure;
corrValfsimsn = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValfsimsnSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValfsimsn,'w')
bar(corrHistVect,corrValfsimsnSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])
set(gca,'TickDir','out');
xlabel('Correlation Coefficient')
ylabel('Number of Pairwise Comparisons')


spikeGraphName = 'allDistTuningCorrelationFSIMSN';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%MSN vs ChAT, close
selDist = 150;
selInterType = 5;


findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
figure
corrValmsnchat = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValmsnchatSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValmsnchat)
bar(corrHistVect,corrValmsnchatSig,'k')
xlim([-1 1])


%plot cumdists?
figure
hold on
plot(corrHistVect,cumsum(corrValmsnmsn)/sum(corrValmsnmsn),'k')
plot(corrHistVect,cumsum(corrValfsifsi)/sum(corrValfsifsi),'c')
plot(corrHistVect,cumsum(corrValfsimsn)/sum(corrValfsimsn),'r')

%try plotting by distance for first 150. 
selInterType = 1;
for i = 1:9
    selDist = (i-1)*25;
    tempTar = intersect(find(distValStore == selDist),find(interTypeStore == selInterType));
    distCorrHistStore(:,i) = hist(corrValStore(tempTar),corrHistVect);
    normDistCorrHistStore(:,i) = distCorrHistStore(:,i)/max(distCorrHistStore(:,i));
end

figure
for i = 1:9
    subplot(9,1,i)
    plot(corrHistVect,normDistCorrHistStore(:,i));
end

%now lets try and dig into the units that I'm pulling out. Lets pull out
%significant MSN-FSI and split between negative and positive. 
selInterType = 1;
selDist = 200;
findSigCorrNeg = find(corrValSigStore == 1 & interTypeStore == selInterType & corrValStore < 0 & distValStore < selDist);
findSigCorrPos = find(corrValSigStore == 1 & interTypeStore == selInterType & corrValStore > 0 & distValStore < selDist);

%now we need to look at cellTrackerStore and extract the MSNs
for i = 1:length(findSigCorrNeg)
    testCells = cellTrackerStore(findSigCorrNeg(i),:);
    if testCells(4) == 0
        sigCorrNegCells(i) = testCells(2);
    elseif testCells(5) == 0
        sigCorrNegCells(i) = testCells(3);
    end
end

for i = 1:length(findSigCorrPos)
    testCells = cellTrackerStore(findSigCorrPos(i),:);
    if testCells(4) == 0
        sigCorrPosCells(i) = testCells(2);
    elseif testCells(5) == 0
        sigCorrPosCells(i) = testCells(3);
    end
end
%this process can generate duplicates if a single cell has multiple
%negative correlations with neighboring FSIs. Eliminate duplicates. 
sigCorrNegCells = unique(sigCorrNegCells);
sigCorrPosCells = unique(sigCorrPosCells);


%now lets plot out tuning curves
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.03 0.03], [0.03 0.03]);
axisSize = ceil(sqrt(length(sigCorrNegCells)));
for i = 1:length(sigCorrNegCells)
    subplot(axisSize,axisSize,i)
    imagesc(squeeze(binValBigStore(:,:,sigCorrNegCells(i)))')
    colormap('parula')
end

hFig = figure;

axisSize = ceil(sqrt(length(sigCorrPosCells)));
for i = 1:length(sigCorrPosCells)
    subplot(axisSize,axisSize,i)
    imagesc(squeeze(binValBigStore(:,:,sigCorrPosCells(i)))')
    colormap('parula')
end


%seem to see some holes in the negative corr, and some positive responses
%in the positive. This is good news I suppose?


%now lets look at the MSN-MSN interactions of the targeted cells. 
selDist = 200;
selType = 0;
negCellmsnCorr = [];
negCellmsnCorrSig = [];
counter = 1;
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
for i = 1:length(sigCorrNegCells)
    findRow1 = find(cellTrackerStore(:,2) == sigCorrNegCells(i));
    findRow2 = find(cellTrackerStore(:,3) == sigCorrNegCells(i));
    
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

counter = 1;
posCellmsnCorr = [];
posCellmsnCorrSig = [];
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
for i = 1:length(sigCorrPosCells)
    findRow1 = find(cellTrackerStore(:,2) == sigCorrPosCells(i));
    findRow2 = find(cellTrackerStore(:,3) == sigCorrPosCells(i));
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

hFig = figure;
hold on
negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,negCellBarPlot)
bar(corrHistVect,negCellSigBarPlot,'k')
xlim([-1 1])

hFig = figure;
hold on
posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,posCellBarPlot)
bar(corrHistVect,posCellSigBarPlot,'k')
xlim([-1 1])

%now lets do the same, but eliminate possibility for double counting
%interactions, which may exist now. 
selDist = 200;
selType = 0;
negCellmsnCorr = [];
negCellmsnCorrSig = [];
counter = 1;
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
tempTrack = cellTrackerStore;
for i = 1:length(sigCorrNegCells)
    findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
    findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

counter = 1;
posCellmsnCorr = [];
posCellmsnCorrSig = [];
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
tempTrack = cellTrackerStore;
for i = 1:length(sigCorrPosCells)
    findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
    findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

hFig = figure;
hold on
negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,negCellBarPlot)
bar(corrHistVect,negCellSigBarPlot,'k')
xlim([-1 1])

hFig = figure;
hold on
posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,posCellBarPlot)
bar(corrHistVect,posCellSigBarPlot,'k')
xlim([-1 1])

%now lets try and do it only for cells within the same group. 
selDist = 200;
negCellmsnCorr = [];
negCellmsnCorrSig = [];
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
tempTrack = cellTrackerStore;
counter = 1;
for i = 1:length(sigCorrNegCells)
    findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
    if findRow1
        otherVals = tempTrack(findRow1,3);
        otherVals = find(ismember(otherVals,sigCorrNegCells));
        findRow1 = findRow1(otherVals);
    end
    
    
    findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
    if findRow2
        otherVals = tempTrack(findRow2,2);
        otherVals = find(ismember(otherVals,sigCorrNegCells));
        findRow2 = findRow2(otherVals);
    end
    
    
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

counter = 1;
posCellmsnCorr = [];
posCellmsnCorrSig = [];
tempTrack = cellTrackerStore;
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
for i = 1:length(sigCorrPosCells)
    findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
    if findRow1
        otherVals = tempTrack(findRow1,3);
        otherVals = find(ismember(otherVals,sigCorrPosCells));
        findRow1 = findRow1(otherVals);
    end
    
    
    findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
    if findRow2
        otherVals = tempTrack(findRow2,2);
        otherVals = find(ismember(otherVals,sigCorrPosCells));
        findRow2 = findRow2(otherVals);
    end
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

hFig = figure;
hold on
negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,negCellBarPlot)
bar(corrHistVect,negCellSigBarPlot,'k')
xlim([-1 1])

hFig = figure;
hold on
posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,posCellBarPlot)
bar(corrHistVect,posCellSigBarPlot,'k')
xlim([-1 1])

%now lets try and do it only for cells not in the same group. 
selDist = 200;
negCellmsnCorr = [];
negCellmsnCorrSig = [];
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
tempTrack = cellTrackerStore;
counter = 1;
for i = 1:length(sigCorrNegCells)
    findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
    if findRow1
        otherVals = tempTrack(findRow1,3);
        otherVals = find(ismember(otherVals,sigCorrPosCells));
        findRow1 = findRow1(otherVals);
    end
    
    
    findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
    if findRow2
        otherVals = tempTrack(findRow2,2);
        otherVals = find(ismember(otherVals,sigCorrPosCells));
        findRow2 = findRow2(otherVals);
    end
    
    
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

counter = 1;
posCellmsnCorr = [];
posCellmsnCorrSig = [];
tempTrack = cellTrackerStore;
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
for i = 1:length(sigCorrPosCells)
    findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
    if findRow1
        otherVals = tempTrack(findRow1,3);
        otherVals = find(ismember(otherVals,sigCorrNegCells));
        findRow1 = findRow1(otherVals);
    end
    
    
    findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
    if findRow2
        otherVals = tempTrack(findRow2,2);
        otherVals = find(ismember(otherVals,sigCorrNegCells));
        findRow2 = findRow2(otherVals);
    end
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

hFig = figure;
hold on
negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,negCellBarPlot)
bar(corrHistVect,negCellSigBarPlot,'k')
xlim([-1 1])

hFig = figure;
hold on
posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,posCellBarPlot)
bar(corrHistVect,posCellSigBarPlot,'k')
xlim([-1 1])

%now lets try and do it only for cells not in either group. 
selDist = 200;
negCellmsnCorr = [];
negCellmsnCorrSig = [];
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
tempTrack = cellTrackerStore;
counter = 1;
for i = 1:length(sigCorrNegCells)
    findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
    if findRow1
        otherVals = tempTrack(findRow1,3);
        otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
        findRow1 = findRow1(otherVals);
    end
    
    
    findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
    if findRow2
        otherVals = tempTrack(findRow2,2);
        otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
        findRow2 = findRow2(otherVals);
    end
    
    
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

counter = 1;
posCellmsnCorr = [];
posCellmsnCorrSig = [];
tempTrack = cellTrackerStore;
findRow3 = find(distValStore < selDist)';
findRow4 = find(interTypeStore == selType)';
findRow3 = intersect(findRow3,findRow4);
for i = 1:length(sigCorrPosCells)
    findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
    if findRow1
        otherVals = tempTrack(findRow1,3);
        otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
        findRow1 = findRow1(otherVals);
    end
    
    
    findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
    if findRow2
        otherVals = tempTrack(findRow2,2);
        otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
        findRow2 = findRow2(otherVals);
    end
%     findRow3 = find(distValStore < selDist)';
    allFinds = sort([findRow1;findRow2]);
    allFinds = intersect(allFinds,findRow3);
    tempTrack(allFinds,:) = 0;
    posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
    posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
    counter = counter + length(allFinds);
end

hFig = figure;
hold on
negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,negCellBarPlot)
bar(corrHistVect,negCellSigBarPlot,'k')
xlim([-1 1])

hFig = figure;
hold on
posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
bar(corrHistVect,posCellBarPlot)
bar(corrHistVect,posCellSigBarPlot,'k')
xlim([-1 1])

%% Lets look at the waveforms a bit more carefully
for i = 1:size(interpWaveStore,1)
    testWave = interpWaveStore(i,:);
    [pks locs] = max(testWave(100:150));
    maxPoint = locs+100-1;
    diffWave = diff(testWave);
    %get half-width
    firstHalf = find(testWave(1:maxPoint) > pks/2,1,'first');
    backHalf = find(testWave(firstHalf:end) < pks/2,1,'first');
    halfWidth(i) = backHalf;
    %get peak trough
    try
        minPoint = find(diffWave(maxPoint:end) > 0,1,'first');
        widthVal(i) = minPoint;
    catch
        disp('No Zero Crossing')
        [pks locs] = min(testWave(maxPoint:end));
        widthVal(i) = locs;
    end
end

figure
stem3(widthVal/300000,halfWidth/300000,bigMaster(:,8))


%% Find cell types
%now plot things!
[indCellType] = functionCellStringFind(masterHeader,'CellType');


%determine if there are units in each category
% findPVs = find(bigMaster(:,indCellType) == 1);
% 
% findMSNs = find(bigMaster(:,indCellType) == 0);

findCHATs = find(bigMaster(:,5) > 0.0005 & bigMaster(:,6) < 1.1);


findPVs = find(bigMaster(:,5) < 0.0004 & bigMaster(:,6) > 1.1);
findMSNs = find(bigMaster(:,5) > 0.0005 & bigMaster(:,6) > 1.1); 
badRPV = find(rpvStore > 0.5);
findBads = intersect(findPVs,badRPV);
findPVs(ismember(findPVs,findBads)) = [];
findBads = intersect(findMSNs,badRPV);
findMSNs(ismember(findMSNs,findBads)) = [];

%% Pull tuning widths
%pull widths
pvWidths = widthStore(:,findPVs,:);
msnWidths = widthStore(:,findMSNs,:);
firstWidthMSN = msnWidths(:,:,2);
firstWidthPV = pvWidths(:,:,2);

%calculate widths with height based system
bigWidthHeightStore = [];
for i = 1:length(bigMaster)
    testValues = squeeze(binValBigStore(:,:,i));
    for j = 1:size(testValues,2)
        [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues(:,j),0.1,3);
        bigWidthHeightStore(j,:,i) = widthVals;
        bigWidthMaxPosStore(j,i) = maxPos;
        bigWidthMaxValStore(j,i) = maxVal;
        bigWidthCutVal(j,i) = cutVal;
    end
end

%now lets just pull out the outward search, since this looks to be far more
%accurate
bigWidthSelWidth = bigWidthHeightStore(:,[1,3],:);

%now lets try and plot out things!
counter = 1;
for i = 1:length(bigMaster)
    if ismember(i,[1,101,201,301,401,501,601])
        hFig = figure;
        counter = 1;
    end
    subplot(10,10,counter)
%     plotData = squeeze(binValBigStore(:,5,i));
    plotData = smooth(squeeze(binValBigStore(:,5,i)),3);
    
    if ismember(i,findPVs)
        plot(plotData,'r','LineWidth',2)
    elseif ismember(i,findMSNs)
        plot(plotData,'k','LineWidth',2)
    else
        plot(plotData,'Color',[0.7 0.7 0.7],'LineWidth',2)
    end
    hold on
    plot([1 16],[0 0],'k')
    plot([1 16],[bigWidthCutVal(5,i) bigWidthCutVal(5,i)],'c','LineWidth',2)
    plot([bigWidthMaxPosStore(5,i) - bigWidthSelWidth(5,1,i) bigWidthMaxPosStore(5,i) - bigWidthSelWidth(5,1,i)],[0 max(squeeze(binValBigStore(:,5,i)))],'g','LineWidth',2)
    plot([bigWidthMaxPosStore(5,i) + bigWidthSelWidth(5,2,i) bigWidthMaxPosStore(5,i) + bigWidthSelWidth(5,2,i)],[0 max(squeeze(binValBigStore(:,5,i)))],'g','LineWidth',2)
    plot([bigWidthMaxPosStore(5,i) bigWidthMaxPosStore(5,i)],[0 max(squeeze(binValBigStore(:,5,i)))],'m','LineWidth',2)
    xlim([1 16])
    counter = counter + 1;
end

%just plot FSIs
hFig = figure;
tarDB = 1;
plotScale = ceil(sqrt(length(findPVs)));
for i = 1:length(findPVs)
    subplot(10,10,i)
%     plotData = squeeze(binValBigStore(:,5,i));
    plotData = smooth(squeeze(binValBigStore(:,tarDB,findPVs(i))),3);
    plot(plotData,'r','LineWidth',2)
    hold on
    plot([1 16],[0 0],'k')
    plot([1 16],[bigWidthCutVal(tarDB,findPVs(i)) bigWidthCutVal(tarDB,findPVs(i))],'c','LineWidth',2)
    plot([bigWidthMaxPosStore(tarDB,findPVs(i)) - bigWidthSelWidth(tarDB,1,findPVs(i)) bigWidthMaxPosStore(tarDB,findPVs(i)) - bigWidthSelWidth(tarDB,1,findPVs(i))],[0 max(squeeze(binValBigStore(:,tarDB,findPVs(i))))],'g','LineWidth',2)
    plot([bigWidthMaxPosStore(tarDB,findPVs(i)) + bigWidthSelWidth(tarDB,2,findPVs(i)) bigWidthMaxPosStore(tarDB,findPVs(i)) + bigWidthSelWidth(tarDB,2,findPVs(i))],[0 max(squeeze(binValBigStore(:,tarDB,findPVs(i))))],'g','LineWidth',2)
    plot([bigWidthMaxPosStore(tarDB,findPVs(i)) bigWidthMaxPosStore(tarDB,findPVs(i))],[0 max(squeeze(binValBigStore(:,tarDB,findPVs(i))))],'m','LineWidth',2)
    xlim([1 16])
end


%% Go through and store specific values. MSNs
%now we need to go through systematically. Store values like threshold
%value, width 10 dB above threshold, shape of response (monotonic
%increasing?)
infoStoreMSN = zeros(length(firstWidthMSN),6);
for i = 1:length(firstWidthMSN)
    %store the total number of DB steps with responses
    infoStoreMSN(i,1) = length(find(firstWidthMSN(:,i) > 0));
    %next, store the "integral" of responses, collapsing across dB
    infoStoreMSN(i,2) = sum(firstWidthMSN(:,i));
    %next, we need to try and look for threshold value. Do this by creeping
    %across. First, check to see if there are any values at all
    if length(find(firstWidthMSN(:,i) > 0)) <=1
        disp('Insufficient Responses Found')
        %report no threshold
        infoStoreMSN(i,3) = 0;
        %report no width
        infoStoreMSN(i,4) = 0;
        %report no shape
        infoStoreMSN(i,5) = 0;
    else
        disp('Sufficient Responses Found')
        %this assumes that there is stuff! now try and find things.
        currVect = firstWidthMSN(:,i);
        %calculate threshold value in different ways, one by the first
        %point after, one at the last zero.
        threshVal1 = find(currVect > 0,1,'first');
        threshVal2 = find(currVect == 0,1,'last')+1;
        if threshVal1 == threshVal2 %if two thresholds are equal, then response is continuous
            threshVal = threshVal1;
            %now check to see if monotonic increasing
            infoStoreMSN(i,3) = threshVal; %store threshold value
            %now store width at level above threshold. If threshold equals
            %lowest value, take that width.
            if threshVal < length(currVect) - 1;
                infoStoreMSN(i,4) = currVect(threshVal+1);
            else
                infoStoreMSN(i,4) = currVect(threshVal);
            end
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStoreMSN(i,5) = 1; %store monotonic increasing shape
            else
                infoStoreMSN(i,5) = 2; %store value as non-monotonic
            end
        elseif threshVal1 == 1 & isempty(threshVal2)
            infoStoreMSN(i,3) = 1;%store threshold as 1
            infoStoreMSN(i,4) = currVect(2); %store width
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStoreMSN(i,5) = 1; %store monotonic increasing shape
            else
                infoStoreMSN(i,5) = 2; %store value as non-monotonic
            end
        else
            threshVal = threshVal2;
            infoStoreMSN(i,3) = threshVal; %store threshold value. This is still useful information
            infoStoreMSN(i,4) = 0;
            infoStoreMSN(i,5) = 2;
        end
    end
end
%% Go through and store specific values. PVs
%do the same for PVs
infoStorePV = zeros(length(firstWidthPV),6);
for i = 1:length(firstWidthPV)
    %store the total number of DB steps with responses
    infoStorePV(i,1) = length(find(firstWidthPV(:,i) > 0));
    %next, store the "integral" of responses, collapsing across dB
    infoStorePV(i,2) = sum(firstWidthPV(:,i));
    %next, we need to try and look for threshold value. Do this by creeping
    %across. First, check to see if there are any values at all
    if length(find(firstWidthPV(:,i) > 0)) <=1
        disp('Insufficient Responses Found')
        %report no threshold
        infoStorePV(i,3) = 0;
        %report no width
        infoStorePV(i,4) = 0;
        %report no shape
        infoStorePV(i,5) = 0;
    else
        disp('Sufficient Responses Found')
        %this assumes that there is stuff! now try and find things.
        currVect = firstWidthPV(:,i);
        %calculate threshold value in different ways, one by the first
        %point after, one at the last zero.
        threshVal1 = find(currVect > 0,1,'first');
        threshVal2 = find(currVect == 0,1,'last')+1;
        if threshVal1 == threshVal2 %if two thresholds are equal, then response is continuous
            threshVal = threshVal1;
            %now check to see if monotonic increasing
            infoStorePV(i,3) = threshVal; %store threshold value
            if threshVal < length(currVect) - 1;
                infoStorePV(i,4) = currVect(threshVal+1);
            else
                infoStorePV(i,4) = currVect(threshVal);
            end
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStorePV(i,5) = 1; %store monotonic increasing shape
            else
                infoStorePV(i,5) = 2; %store value as non-monotonic
            end
        elseif threshVal1 == 1 & isempty(threshVal2)
            infoStorePV(i,3) = 1;%store threshold as 1
            infoStorePV(i,4) = currVect(2); %store width
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStorePV(i,5) = 1; %store monotonic increasing shape
            else
                infoStorePV(i,5) = 2; %store value as non-monotonic
            end
        else
            threshVal = threshVal2;
            infoStorePV(i,3) = threshVal; %store threshold value. This is still useful information
            infoStorePV(i,4) = 0;
            infoStorePV(i,5) = 2;
        end
    end
end
%% SAVE COMBINED DATA
condensedMSN = infoStoreMSN(infoStoreMSN(:,5) == 1,:);
condensedPV = infoStorePV(infoStorePV(:,5) == 1,:);

[indPkTr] = functionCellStringFind(masterHeader,'PeakTrough');
[indISI] = functionCellStringFind(masterHeader,'isiCov');
[indBaseFire] = functionCellStringFind(masterHeader,'BaselineRate');
[indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
[indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');

%make vector for width. if neuron responds to everything, that makes 5 *
%16, or 80.
widthHistVect = [0:1:40];

%now lets generate half-width values

for i = 1:length(bigMaster)
    testWave = interpWaveStore(i,:);
    %determine max
    maxVal = max(testWave(1:150));
    halfFirst = find(testWave(1:150) > maxVal/2,1,'first');
    halfSecond = find(testWave(halfFirst:end) > maxVal/2,1,'last');
    halfWidth(i) = halfSecond/300000;
end

%% first figure, OVERVIEW OF FIRING RATES AND PEAK TROUGH, RESPONSE TYPES ETC
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.1 0.1], [0.1 0.1]);
hFig = figure;
set(hFig, 'Position', [10 80 500 1000])

%column 1 plot spike width vs coefficient of variation
subplot(3,1,1)
hold on

plot(bigMaster(:,indPkTr),bigMaster(:,indISI),'.','Color',[0.7 0.7 0.7])
plot(bigMaster(bigMaster(:,indCellType) == 0,indPkTr),bigMaster(bigMaster(:,indCellType) == 0,indISI),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),bigMaster(bigMaster(:,indCellType) == 1,indISI),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),bigMaster(bigMaster(:,indCellType) == 2,indISI),'g.')
xlabel('Peak Trough (ms)')
ylabel('ISI Coefficient of Variation')
set(gca,'TickDir','out');

%now plot out spike width vs half width
subplot(3,1,2)
hold on

plot(bigMaster(:,indPkTr),halfWidth(:),'.','Color',[0.7 0.7 0.7])
plot(bigMaster(bigMaster(:,indCellType) == 0,indPkTr),halfWidth(bigMaster(:,indCellType) == 0),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),halfWidth(bigMaster(:,indCellType) == 1),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),halfWidth(bigMaster(:,indCellType) == 2),'g.')
xlabel('Peak Trough (ms)')
ylabel('Half-Width (ms)')
set(gca,'TickDir','out');

%plot out spike width with FR
subplot(3,1,3)
hold on

plot(bigMaster(:,indPkTr),bigMaster(:,indBaseFire),'.','Color',[0.7 0.7 0.7])
plot(bigMaster(bigMaster(:,indCellType) == 0,indPkTr),bigMaster(bigMaster(:,indCellType) == 0,indBaseFire),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),bigMaster(bigMaster(:,indCellType) == 1,indBaseFire),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),bigMaster(bigMaster(:,indCellType) == 2,indBaseFire),'g.')
xlabel('Peak Trough (ms)')
ylabel('Baseline Firing Rate')
set(gca,'TickDir','out');

spikeGraphName = 'WaveformAnalysis';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')




subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1000 1000])
%plot out pie chart of difference cells

cellDist = [length(findMSNs),length(findPVs),length(findCHATs)];
pie(cellDist)
labels = {strcat('MSNs(',num2str(length(findMSNs)),')'),strcat('PVs(',num2str(length(findPVs)),')'),strcat('ChATs(',num2str(length(findCHATs)),')')};
detZero = find(cellDist == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')

spikeGraphName = 'CellTypePieChart';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%plot out MSN and FSI pie charts
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])

subplot(2,2,1)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findMSNs,:);
plot(waveHolder')
title('MSN Waveforms')

%plot distribution of response types
subplot(2,2,3)
holder = bigMaster(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('MSNs n=',num2str(length(findMSNs))))

%plot out PV stuff

subplot(2,2,2)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findPVs,:);
plot(waveHolder')
title('PV Waveforms')

%plot distribution of response types
subplot(2,2,4)
holder = bigMaster(findPVs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('PVs n=',num2str(length(findPVs))))

spikeGraphName = 'PVMSNPieCharts';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')




subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])

%column 1 plot spike width vs coefficient of variation
subplot(4,4,1)
hold on

plot(bigMaster(:,indPkTr),bigMaster(:,indISI),'.','Color',[0.7 0.7 0.7])
plot(bigMaster(bigMaster(:,indCellType) == 0,indPkTr),bigMaster(bigMaster(:,indCellType) == 0,indISI),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),bigMaster(bigMaster(:,indCellType) == 1,indISI),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),bigMaster(bigMaster(:,indCellType) == 2,indISI),'g.')
xlabel('Peak Trough (ms)')
ylabel('ISI Coefficient of Variation')

%now plot out spike width vs half width
subplot(4,4,5)
hold on

plot(bigMaster(:,indPkTr),halfWidth(:),'.','Color',[0.7 0.7 0.7])
plot(bigMaster(bigMaster(:,indCellType) == 0,indPkTr),halfWidth(bigMaster(:,indCellType) == 0),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),halfWidth(bigMaster(:,indCellType) == 1),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),halfWidth(bigMaster(:,indCellType) == 2),'g.')
xlabel('Peak Trough (ms)')
ylabel('Half-Width (ms)')

%plot out spike width with FR
subplot(4,4,9)
hold on

plot(bigMaster(:,indPkTr),bigMaster(:,indBaseFire),'.','Color',[0.7 0.7 0.7])
plot(bigMaster(bigMaster(:,indCellType) == 0,indPkTr),bigMaster(bigMaster(:,indCellType) == 0,indBaseFire),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),bigMaster(bigMaster(:,indCellType) == 1,indBaseFire),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),bigMaster(bigMaster(:,indCellType) == 2,indBaseFire),'g.')
xlabel('Peak Trough (ms)')
ylabel('Baseline Firing Rate')

%plot out pie chart of difference cells
subplot(4,4,13)
cellDist = [length(findMSNs),length(findPVs),length(findCHATs)];
pie(cellDist)
labels = {strcat('MSNs(',num2str(length(findMSNs)),')'),strcat('PVs(',num2str(length(findPVs)),')'),strcat('ChATs(',num2str(length(findCHATs)),')')};
detZero = find(cellDist == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')

%column 2, plot out pie charts of MSNs

subplot(4,4,2)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findMSNs,:);
plot(waveHolder')
title('MSN Waveforms')

%plot distribution of response types
subplot(4,4,6)
holder = bigMaster(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('MSNs n=',num2str(length(findMSNs))))

%plot out integrals of response plots
subplot(4,4,10)
holder = intFastPos(findMSNs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))


subplot(4,4,14)
hold on
holder = intSlowPos(findMSNs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))


%column 3, plot out PV stuff

subplot(4,4,3)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findPVs,:);
plot(waveHolder')
title('PV Waveforms')

%plot distribution of response types
subplot(4,4,7)
holder = bigMaster(findPVs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('PVs n=',num2str(length(findPVs))))

%plot out integrals of response plots
subplot(4,4,11)
holder = intFastPos(findPVs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

subplot(4,4,15)
hold on
holder = intSlowPos(findPVs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

%column 4, plot out pie charts of CHATs

subplot(4,4,4)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findCHATs,:);
plot(waveHolder')
title('CHAT Waveforms')

%plot distribution of response types
subplot(4,4,8)
holder = bigMaster(findCHATs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('CHATs n=',num2str(length(findCHATs))))

%plot out integrals of response plots
subplot(4,4,12)
holder = intFastPos(findCHATs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

subplot(4,4,16)
hold on
holder = intSlowPos(findCHATs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))


spikeGraphName = 'WrapperFigure1Overview';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Generate tuning curves that are step functions based on significance.
%how lets look at sigValBigStore
sigValConv = sigValBigStore;
sigValConv(sigValConv <= 0.001) = 4;
sigValConv(sigValConv <= 0.01) = 3;
sigValConv(sigValConv <= 0.05) = 2;
sigValConv(sigValConv <= 1) = 1;
clims = [1 3];

%% Plot out tuning curves without white noise, actual values and step function significance
%just do PVs
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
plotWid = ceil(sqrt(length(findPVs)));
for i = 1:length(findPVs)
    subplot(plotWid,plotWid,i)
    imagesc(binValBigStore(:,:,findPVs(i))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'pvBinStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findPVs)
    subplot(plotWid,plotWid,i)
    imagesc(sigValConv(:,:,findPVs(i))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end
spikeGraphName = 'pvSigStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%and just MSNs
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
plotWid = ceil(sqrt(length(findMSNs)));
for i = 1:length(findMSNs)
    subplot(plotWid,plotWid,i)
    imagesc(binValBigStore(:,:,findMSNs(i))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'msnBinStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findMSNs)
    subplot(plotWid,plotWid,i)
    imagesc(sigValConv(:,:,findMSNs(i))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'msnSigStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Now lets select for units with five significant positive responses. 
signStore = sign(binValBigStore);
signSigStore = sigValBigStore.*signStore;
signSigStore(signSigStore <= 0) = NaN;
sigCut = 0.01;
for i = 1:size(signSigStore,3)
    findSigs(i) = length(find(signSigStore(:,:,i) <= sigCut));
end

figure
hist(findSigs,[0:1:20])
xlim([-0.5 20.5])

%target units with greater than 5 significant responses.
tarCells = find(findSigs >5);


%first, we want to eliminate non-significant latencies
sigMin = 0.01; %minimum significance value
binSig = sigValBigStore;
binSig(binSig > sigMin) = NaN;
binSig(binSig <= sigMin) = 1;

latConv = binSig.*latMapBigStore;
widthLatConv = binSig.*widthLatStore;
%now delete zero latency values
latConv(latConv == 0) = NaN;

%lets separate white noise from the rest
% latConvWhite = latConv(1,:,:);
latConvTone = latConv(:,:,:);
latConvWidthTone = widthLatConv(:,:,:);


%now lets extract targeted latencies. THIS IS FOR BIGGEST RESPONSE
for i = 1:length(tarCells);
    tempLat = latConvTone(:,:,tarCells(i));
    tarLats(i) = tempLat(bigFreqStore(tarCells(i)),bigDBStore(tarCells(i)));
end
%QC check, remove all latencies and tarCells that are NaNs
nanFind = isnan(tarLats);
tarLats(nanFind) = [];
tarCells(nanFind) = [];

[C,tarPVs,ib] = intersect(tarCells,findPVs);
[C,tarMSNs,ib] = intersect(tarCells,findMSNs);

latHistVect = [0:0.001:0.1];

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(3,2,1)
hist(tarLats(tarPVs),latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(tarPVs)), ' FSI Min Latency Tone'))
set(gca,'TickDir','out');

% subplot(3,2,2) hist(minLatTonePV,latHistVect) xlim([latHistVect(1)
% latHistVect(end)])
% title(strcat(num2str(length(find(~isnan(minLatTonePV)))),' FSI Min
% Latency Pure Tone')) % title('FSI Min Latency Pure Tone')

subplot(3,2,3)
hist(tarLats(tarMSNs),latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(tarMSNs)),' MSN Min Latency Tone'))
set(gca,'TickDir','out');
% title('MSN Min Latency White Noise')

% subplot(3,2,4) hist(minLatToneMSN,latHistVect) xlim([latHistVect(1)
% latHistVect(end)])
% title(strcat(num2str(length(find(~isnan(minLatToneMSN)))),' MSN Min
% Latency Pure Tone')) % title('MSN Min Latency Pure Tone')

subplot(3,1,3)
hold on
bar(1:2,[nanmean(tarLats(tarPVs)),nanmean(tarLats(tarMSNs))],'w')
errorbar(1:2,[nanmean(tarLats(tarPVs)),nanmean(tarLats(tarMSNs))],[nanstd(tarLats(tarPVs)),nanstd(tarLats(tarMSNs))])
set(gca,'TickDir','out');
% xticks([1:4]) xticklabels({'FSI White','MSN White','FSI Tone','MSN
% Tone'})

spikeGraphName = '5ValBigRespLat';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at just 70dB
% 
% minLatWhitePV = squeeze((latConvWhite(:,5,tarCells(tarPVs))));
% minLatWhiteMSN = squeeze((latConvWhite(:,5,tarCells(tarMSNs))));
minLatTonePV = squeeze((min(latConvTone(:,5,tarCells(tarPVs)))));
minLatToneMSN = squeeze((min(latConvTone(:,5,tarCells(tarMSNs)))));

latHistVect = [0:0.001:0.1];

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
% subplot(3,2,1)
% hist(minLatWhitePV,latHistVect)
% xlim([latHistVect(1) latHistVect(end)])
% title(strcat(num2str(length(find(~isnan(tarPVs)))), ' FSI Min Latency White Noise'))
% 
subplot(3,2,2)
hist(minLatTonePV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(tarPVs)))),' FSI Min Latency Pure Tone'))
set(gca,'TickDir','out');
% title('FSI Min Latency Pure Tone')

% subplot(3,2,3)
% hist(minLatWhiteMSN,latHistVect)
% xlim([latHistVect(1) latHistVect(end)])
% title(strcat(num2str(length(find(~isnan(tarMSNs)))),' MSN Min Latency White Noise'))
% % title('MSN Min Latency White Noise')

subplot(3,2,4)
hist(minLatToneMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(tarMSNs)))),' MSN Min Latency Pure Tone'))
set(gca,'TickDir','out');
% title('MSN Min Latency Pure Tone')

subplot(3,1,3)
hold on
bar(1:2,[nanmean(minLatTonePV),nanmean(minLatToneMSN)],'w')
errorbar(1:2,[nanmean(minLatTonePV),nanmean(minLatToneMSN)],[nanstd(minLatTonePV),nanstd(minLatToneMSN)])
set(gca,'TickDir','out');
% xticks([1:4]) xticklabels({'FSI White','MSN White','FSI Tone','MSN
% Tone'})

spikeGraphName = '5ValMinLat';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%the results of this demonstrate that there is a weak difference in latency
%that is insignificant. Looking at histograms, looks like this might be in
%part due to a number of slower responding PV cells.

%lets try and look on a recording by recording basis
tarRec = recStore(tarCells)';
tarRec(:,2) = NaN;
tarRec(tarPVs,2) = 1;
tarRec(tarMSNs,2) = 0;
tarRec(:,3) = 0;
%check to see if specific recording has both PVs and MSNs
for i = 1:numFiles
    tempStore = tarRec(tarRec(:,1) == i,2);
    if sum(tempStore) == 0 | sum(tempStore) == length(tempStore)
        disp('Only One Cell Type')
    else
        disp('Multiple Cell Types')
        tarRec(tarRec(:,1) == i,3) = i;
    end
end

tarRec(:,4) = tarLats;
tester = [-.2:0.0005:0.4];
%Doesnt look like a per-recording analysis will pull out anything different
%really.



%now lets plot out heatmaps?

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
plotWid = ceil(sqrt(length(tarPVs)));
for i = 1:length(tarPVs)
    subplot(plotWid,plotWid,i)
    imagesc(binValBigStore(:,:,tarCells(tarPVs(i)))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = '5ValFSIBinStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
plotWid = ceil(sqrt(length(tarMSNs)));
for i = 1:length(tarMSNs)
    subplot(plotWid,plotWid,i)
    imagesc(binValBigStore(:,:,tarCells(tarMSNs(i)))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = '5ValMSNBinStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets z score everything first
tarHists = fineHist(tarCells,:);
for i = 1:length(tarCells)
    %extract std and mean from baseline
    meanVal = mean(tarHists(i,1:401));
    stdVal = std(tarHists(i,1:401));
    zHists(i,:) = (tarHists(i,:) - meanVal)/stdVal;
end

hFig = figure;
hold on
plot([-.2:0.0005:0.4],mean(zHists((tarPVs),:))-mean(zHists((tarPVs),401)),'r')
plot([-.2:0.0005:0.4],mean(zHists((tarMSNs),:))-mean(zHists((tarMSNs),401)),'k')
set(gca,'TickDir','out');
spikeGraphName = '5ValAverageZPlotsMSNrFSIb';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
hold on
plot([-.2:0.0005:0.4],mean(zHists((tarPVs),:))-mean(zHists((tarPVs),401)),'r')
plot([-.2:0.0005:0.4],mean(zHists((tarMSNs),:))-mean(zHists((tarMSNs),401)),'k')
set(gca,'TickDir','out');
% xlim([-0.02 0.05])
xlim([0 0.02])
spikeGraphName = '5ValAverageZPlotsMSNrFSIbZOOM';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
hold on
for i = 1:length(tarPVs)
    plot([-.2:0.0005:0.4],zHists((tarPVs(i)),:),'r')
end
for i = 1:length(tarMSNs)
    plot([-.2:0.0005:0.4],zHists((tarMSNs(i)),:),'k')
end
xlim([0 0.02])
set(gca,'TickDir','out');
spikeGraphName = '5ValPVMSNzOverallPlotIndiv';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%we can also now try and do this for all units. lets z score everything
%first
% tarHists = fineHist;
latFineHist(:,1) = 0;
latFineHist(:,end) = 0;
for i = 1:length(tarCells)
    %extract std and mean from baseline
    meanVal = mean(latFineHist(tarCells(i),1:401));
    stdVal = std(latFineHist(tarCells(i),1:401));
    allZHists(i,:) = (latFineHist(tarCells(i),:) - meanVal)/stdVal;
end

%this has been fixed now. however, still some shitty units. lets try and
%evaluate which ones can be salvaged.
for i = 1:length(tarCells)
    %extract std and mean from baseline
    toneSpikes(i) = sum(latFineHist(tarCells(i),401:601))-sum(latFineHist(tarCells(i),1:401));
    peakVals(i) = max(latFineHist(tarCells(i),401:601))-mean(latFineHist(tarCells(i),1:401));
    [pks maxInd(i)] = max(latFineHist(tarCells(i),401:601));
end

specFind = find(toneSpikes(tarPVs)>0);
hFig = figure;
hold on
for i = 1:length(specFind)
%     plot([-.2:0.0005:0.4],smooth(allZHists((tarPVs(specFind(i))),:),11)/max(smooth(allZHists((tarPVs(specFind(i))),:),11))+i,'Color',[rand(1),rand(1),rand(1)])
    plot(smooth(allZHists((tarPVs(specFind(i))),1:maxInd(tarPVs(specFind(i)))+400),11)/max(smooth(allZHists((tarPVs(specFind(i))),1:maxInd(tarPVs(specFind(i)))+400),11))+i,'Color',[rand(1),rand(1),rand(1)])
end

specFind = find(toneSpikes(tarMSNs) > 10);
hFig = figure;
hold on
for i = 1:length(specFind)
%     plot([-.2:0.0005:0.4],smooth(allZHists((tarMSNs(specFind(i))),:),11)/max(smooth(allZHists((tarMSNs(specFind(i))),:),11))+i,'Color',[rand(1),rand(1),rand(1)])
    plot(smooth(allZHists((tarMSNs(specFind(i))),1:maxInd(tarMSNs(specFind(i)))+400),11)/max(smooth(allZHists((tarMSNs(specFind(i))),1:maxInd(tarMSNs(specFind(i)))+400),11))+i,'Color',[rand(1),rand(1),rand(1)])
end

%lets try and do this without spreading in y axis

specFind = find(toneSpikes(tarPVs)>0);
hFig = figure;
hold on
for i = 1:length(specFind)
%     plot([-.2:0.0005:0.4],smooth(allZHists((tarPVs(specFind(i))),:),11)/max(smooth(allZHists((tarPVs(specFind(i))),:),11))+i,'Color',[rand(1),rand(1),rand(1)])
    plot([-0.2:0.0005:-0.2+0.0005*(maxInd(tarPVs(specFind(i)))+400-1)],smooth(allZHists((tarPVs(specFind(i))),1:maxInd(tarPVs(specFind(i)))+400),11)/max(smooth(allZHists((tarPVs(specFind(i))),1:maxInd(tarPVs(specFind(i)))+400),11)),'r')
end

specFind = find(toneSpikes(tarMSNs) > 10);
% hFig = figure;
hold on
for i = 1:length(specFind)
%     plot([-.2:0.0005:0.4],smooth(allZHists((tarMSNs(specFind(i))),:),11)/max(smooth(allZHists((tarMSNs(specFind(i))),:),11))+i,'Color',[rand(1),rand(1),rand(1)])
    plot([-0.2:0.0005:-0.2+0.0005*(maxInd(tarMSNs(specFind(i)))+400-1)],smooth(allZHists((tarMSNs(specFind(i))),1:maxInd(tarMSNs(specFind(i)))+400),11)/max(smooth(allZHists((tarMSNs(specFind(i))),1:maxInd(tarMSNs(specFind(i)))+400),11)),'k')
end

xlim([0 0.05])
set(gca,'TickDir','out');
spikeGraphName = '5ValPVMSNzSelFreqs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and plot things out so that we get the calculated latency
%from these selected frequencies. Since we've already z-scored things, its
%about selecting a cutoff and using that!

zcut = 2;
smoothZ = [];
for i = 1:size(allZHists,1);
    smoothZ(i,:) = smooth(allZHists(i,:),11);
end

for i = 1:size(smoothZ,1);
    testData = smoothZ(i,:);
    try
        findFirst = find(testData(400:end) >= zcut,1,'first');
        latStore(i) = findFirst * binSize;
    catch
        latStore(i) = NaN;
    end
end


hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(2,1,1)
hist(latStore(tarMSNs),latHistVect)
xlim([latHistVect(1) latHistVect(end)])
set(gca,'TickDir','out');
title('MSN Latency Selected Freqs 5ms Smooth')

subplot(2,1,2)
hist(latStore(tarPVs),latHistVect)
xlim([latHistVect(1) latHistVect(end)])
set(gca,'TickDir','out');
title('FSI Latency Selected Freqs 5ms Smooth')
% title('FSI Min Latency Pure Tone')

spikeGraphName = '5ValPVMSNLatHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%now lets pull the height based width values, and the centers of these
%"responses"
% pvHeightWidth = bigWidthSelWidth(:,:,tarCells(tarPVs));
% msnHeightWidth = bigWidthSelWidth(:,:,tarCells(tarMSNs));
tarCellHeightWidth = bigWidthSelWidth(:,:,tarCells);

% pvMaxPos = bigWidthMaxPosStore(:,tarCells(tarPVs));
% msnMaxPos = bigWidthMaxPosStore(:,tarCells(tarMSNs));
tarCellMaxPos = bigWidthMaxPosStore(:,tarCells);

%now lets determine whether significant values exist at each DB level.
sigCut = 0.05;
abbrevSig = signSigStore(:,:,tarCells);
for i = 1:length(tarCells)
    tester = min(squeeze(abbrevSig(:,:,i)));
    tester(tester <= sigCut) = 10;
    tester(tester < 10) = 0;
    tester(tester == 10) = 1;
    
    sigStore5Val(:,i) = tester;
end
sigStore5Val(sigStore5Val~=1) = 0;

%ok, now we have significant values, max points, and the detected width in
%each direction. Lets try and pull things out now?

sigWidthStore = zeros(size(sigStore5Val));
sideWarnLow = zeros(size(sigStore5Val));
sideWarnHi = zeros(size(sigStore5Val));
for i = 1:length(tarCells)
    for j = 1:size(sigStore5Val,1)
        %first, check if there is a significant value
        if sigStore5Val(j,i) == 1
            %now check for width values. If NaN on either side, determine
            %edge.
            %start with lower edge. 
            if isnan(tarCellHeightWidth(j,1,i))
                lowWidth = tarCellMaxPos(j,i);
                sideWarnLow(j,i) = 1;
            else
                lowWidth = tarCellHeightWidth(j,1,i);
            end
            if isnan(tarCellHeightWidth(j,2,i))
                hiWidth = 16 - tarCellMaxPos(j,i);
                sideWarnHi(j,i) = 1;
            else
                hiWidth = tarCellHeightWidth(j,2,i);
            end
            sigWidthStore(j,i) = hiWidth + lowWidth;
        end
    end
end

%now separate to MSNs vs FSIs
sigWidthPV = sigWidthStore(:,tarPVs);
sigWidthMSN = sigWidthStore(:,tarMSNs);

%now lets just do some ugly stuff and look only at 70dB point. 
widthValsPV = sigWidthPV(5,:);
widthValsMSN = sigWidthMSN(5,:);

widthHistVect = [0:0.5:16];
histValPV = hist(widthValsPV,widthHistVect);
histValMSN = hist(widthValsMSN,widthHistVect);


% now lets try and plot out width, threshold, with the width of 3 setup.
%first do PVs.
threshStorePV = [];
for i = 1:length(tarPVs)
    tester = sigValConv(2:end,:,tarCells(tarPVs(i)));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    threshStorePV(i) = testFind;
end

%now do MSNs
threshStoreMSN = [];
for i = 1:length(tarMSNs)
    tester = sigValConv(2:end,:,tarCells(tarMSNs(i)));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    threshStoreMSN(i) = testFind;
end
%to access width data, need to get correct indices

[C pvIndex b] = intersect(findPVs,tarCells(tarPVs));
[C msnsIndex b] = intersect(findMSNs,tarCells(tarMSNs));

%set ylimits for plots

ylimThresh = [0 35];
ylimWidth = [0 15];
%plot this out.
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,2,1)
hist(threshStorePV,[1:1:5],'r')
xlim([0.5 5.5])
ylim(ylimThresh)
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});
title(strcat('Amplitude Threshold (FSIs) Mean:',num2str(mean(threshStorePV)),'pval',num2str(ranksum(threshStorePV,threshStoreMSN))))
set(gca,'TickDir','out');

subplot(2,2,3)
hist(threshStoreMSN,[1:1:5],'k')
xlim([0.5 5.5])
title(strcat('Amplitude Threshold (MSNs) Mean:',num2str(mean(threshStoreMSN))))
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});
ylim(ylimThresh)
set(gca,'TickDir','out');

subplot(2,2,2)
hist(widthValsPV,[1:1:16],'r')
xlim([0.5 16.5])
ylim(ylimWidth)
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('70 dB Tuning Width (FSIs) Mean:',num2str(mean(widthValsPV)),'pval',num2str(ranksum(widthValsPV,widthValsMSN))))
% title('Tuning Width (FSIs)')
set(gca,'TickDir','out');

subplot(2,2,4)
hist(widthValsMSN,[1:1:16],'k')
xlim([0.5 16.5])
ylim(ylimWidth)
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('70 dB Tuning Width (MSNs) Mean:',num2str(mean(widthValsMSN))))
set(gca,'TickDir','out');

spikeGraphName = '5ValIntensityThreshAndWidth';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

% plot out BFs of width 3 selected units.

hFig = figure;
xlims = [4000 32000];
ylims = [0 15];
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hist(bfStore(tarCells(tarPVs)),[4000;4594.79342000000;5278.03164300000;6062.86626600000;6964.40450600000;8000;9189.58684000000;10556.0632900000;12125.7325300000;13928.8090100000;16000;18379.1736800000;21112.1265700000;24251.4650600000;27857.6180300000;32000],'r')
set(gca,'xscale','log')
xlim(xlims)
ylim(ylims)
% xlim([4000 32000])
ylabel('Number of Cells')
xlabel('Best Frequency (Hz)')
% set(gca,'XTickLabel',{'<=20','30','40','60','70'});
title('BFs (FSIs)')
set(gca,'TickDir','out');

subplot(2,1,2)
hist(bfStore(tarCells(tarMSNs)),[4000;4594.79342000000;5278.03164300000;6062.86626600000;6964.40450600000;8000;9189.58684000000;10556.0632900000;12125.7325300000;13928.8090100000;16000;18379.1736800000;21112.1265700000;24251.4650600000;27857.6180300000;32000],'k')
set(gca,'xscale','log')
xlim(xlims)
ylim(ylims)
% xlim([4000 32000])
ylabel('Number of Cells')
xlabel('Best Frequency (Hz)')
% set(gca,'XTickLabel',{'<=20','30','40','60','70'});
title('BFs (MSNs)')
set(gca,'TickDir','out');


spikeGraphName = '5ValBFPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets plot out width by threshold amplitude.

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [80 80 1900 1000])
for i = 1:5
    subplot(2,5,i)
    %find PVs with threshold of a target value
    finderPV = find(threshStorePV == i);
    finderMSN = find(threshStoreMSN == i);
    if finderPV
        hist(sigWidthPV(i,finderPV),[1:1:16],'r')
        testVal = mean(sigWidthPV(i,(finderPV)));
        title(strcat('FSI Width for Amp Level',num2str(i),'Mean',num2str(testVal),'pval',num2str(ranksum(sigWidthPV(i,finderPV),sigWidthMSN(i,finderMSN)))))
        set(gca,'TickDir','out');
    end
    subplot(2,5,i+5)
    %find MSNs with threshold of a target value
    if finderMSN
        hist(sigWidthMSN(i,(finderMSN)),[1:1:16],'k')
        testVal = mean(sigWidthMSN(i,(finderMSN)));
        title(strcat('MSN Width for Amp Level',num2str(i),'Mean',num2str(testVal)))
        set(gca,'TickDir','out');
    end
end
set(gca,'TickDir','out');

spikeGraphName = '5ValPlotWidthBasedOnThreshAmp';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%now lets isolate the units that are responsive at 20 dB and see how widths
%change over dB


hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [80 80 1900 1000])
targetThresh = 1;
finderPV = find(threshStorePV == targetThresh);
finderMSN = find(threshStoreMSN == targetThresh);
ylims = [0 13];
xlims = [0 16];

for i = 1:5
    subplot(2,5,i)
    %find PVs with threshold of a target value
    if finderPV
        hist(sigWidthPV(i,finderPV),[1:1:16],'r')
        testVal = mean(sigWidthPV(i,(finderPV)));
        title(strcat('FSI Width for Amp Level',num2str(i),'Mean',num2str(testVal),'pval',num2str(ranksum(sigWidthPV(i,finderPV),sigWidthMSN(i,finderMSN)))))
        set(gca,'TickDir','out');
        ylim(ylims)
        xlim(xlims)
    end
    subplot(2,5,i+5)
    %find MSNs with threshold of a target value
    if finderMSN
        hist(sigWidthMSN(i,(finderMSN)),[1:1:16],'k')
        testVal = mean(sigWidthMSN(i,(finderMSN)));
        title(strcat('MSN Width for Amp Level',num2str(i),'Mean',num2str(testVal)))
        set(gca,'TickDir','out');
        ylim(ylims)
        xlim(xlims)
    end
end
set(gca,'TickDir','out');

spikeGraphName = '5ValPlotWidthJust20ThreshAcrossDB';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets plot out individual changes in width for MSNs and PVs
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [80 80 1900 1000])
targetThresh = 1;
finderPV = find(threshStorePV == targetThresh);
finderMSN = find(threshStoreMSN == targetThresh);
hold on
% for i = 1:length(finderPV)
%     plot(sigWidthPV(:,finderPV(i)),'r.-')
% end
errorbar(mean(sigWidthPV(:,finderPV)'),std(sigWidthPV(:,finderPV)')/sqrt(length(finderPV)),'r','LineWidth',2);
% for i = 1:length(finderMSN)
%     plot(sigWidthMSN(:,finderMSN(i)),'k.-')
% end
errorbar(mean(sigWidthMSN(:,finderMSN)'),std(sigWidthMSN(:,finderMSN)')/sqrt(length(finderMSN)),'k','LineWidth',2);
set(gca,'TickDir','out');

spikeGraphName = '5ValPlotMeanWidthChangeAcrossDBs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%now lets plot BF vs width.
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [80 80 1900 1000])
subplot(2,1,1)
plot(bfStore(tarCells(tarPVs)),sigWidthPV(5,:),'k.')
xlabel('BF')
ylabel('Tuning Width')
title('FSI Width vs BF')
set(gca,'TickDir','out');
subplot(2,1,2)
plot(bfStore(tarCells(tarMSNs)),sigWidthMSN(5,:),'k.')
xlabel('BF')
ylabel('Tuning Width')
title('MSN Width vs BF')
set(gca,'TickDir','out');

spikeGraphName = '5ValPlotBFvsWidth';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
 










