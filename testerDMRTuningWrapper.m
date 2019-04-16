
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
        latFineHist(bigMasterInd + j - 1,:) = hist(tarSpikes,[histLims(1):binSize:histLims(2)])/binSize/length(tarStore);

        disp('About to do fine overall histogram')
        %generate a finer scale histogram across all tones
        tempFineRast = functionBasicRaster(s.(desigName{j}).SpikeTimes,s.TrialMatrix(toneFinder,1),histLims);
        fineHistVect = [histLims(1):binSize:histLims(2)];
        fineHist(bigMasterInd + j - 1,:) = hist(tempFineRast(:,1),[histLims(1):binSize:histLims(2)])/binSize/length(toneFinder);
        nameStore{bigMasterInd + j -1} = targetFiles{i};
        unitStore{bigMasterInd + j -1} = desigName{j};
        dmrSpikeNum(bigMasterInd + j - 1) = length(find(s.(desigName{j}).SpikeTimes > s.SoundData.DMRPulses(1) & s.(desigName{j}).SpikeTimes < s.SoundData.DMRPulses(end)));
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
    
    %lets add some code to extract out STA data!
    bigSTAstore(bigMasterInd:bigMasterInd + numUnits - 1,:) = s.STAs;
    bigSTASigstore(bigMasterInd:bigMasterInd + numUnits - 1,:) = s.STASig;
    
    %now lets pull significance information from STA
    newTarget = targetFiles{i};
    newTarget = newTarget(1:end-22);
    newTarget = [newTarget,'DMRData.mat'];
    load(newTarget)
    
    bigStoreSTACorr(bigMasterInd:bigMasterInd + numUnits - 1) = realCorrStore(1:numUnits);
    bigStoreSTACorrBase(bigMasterInd:bigMasterInd + numUnits - 1,:) = corrCoefStore(1:numUnits,:);
    bigStoreSpikeNum(bigMasterInd:bigMasterInd + numUnits - 1) = sum(spikeArray');
    bigMasterInd = bigMasterInd + numUnits;
        
end

faxis = s.DMRfaxis;
toneFreqs = s.SoundData.UniqueFrequencies;

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
figure
corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValmsnmsnSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValmsnmsn)
bar(corrHistVect,corrValmsnmsnSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])
% corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);

%just FSIs, close
selDist = 150;
selInterType = 2;


findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
figure
corrValfsifsi = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValfsifsiSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValfsifsi)
bar(corrHistVect,corrValfsifsiSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])


%MSN vs FSI, close
selDist = 150;
selInterType = 1;


findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
figure
corrValfsimsn = hist(corrValStore(findTarCorr),corrHistVect);
sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
corrValfsimsnSig = hist(corrValStore(sigInt),corrHistVect);
hold on
bar(corrHistVect,corrValfsimsn)
bar(corrHistVect,corrValfsimsnSig,'k')
% hist(corrValStore(findTarCorr),corrHistVect)
xlim([-1 1])

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

%now lets plot out STAs
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.03 0.03], [0.03 0.03]);
axisSize = ceil(sqrt(length(sigCorrNegCells)));
for i = 1:length(sigCorrNegCells)
    subplot(axisSize,axisSize,i)
    imagesc(reshape(bigSTASigstore(sigCorrNegCells(i),:),length(faxis),[]))
    colormap('parula')
end

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.03 0.03], [0.03 0.03]);
axisSize = ceil(sqrt(length(sigCorrPosCells)));
for i = 1:length(sigCorrPosCells)
    subplot(axisSize,axisSize,i)
    imagesc(reshape(bigSTASigstore(sigCorrPosCells(i),:),length(faxis),[]))
    colormap('parula')
end


%plot out overall histograms?
hFig = figure;
hold on
for i = 1:length(sigCorrNegCells)
    plot(smooth(fineHist(sigCorrNegCells(i),:),11)/max(smooth(fineHist(sigCorrNegCells(i),:),11)) + i/2)
end

hFig = figure;
hold on
for i = 1:length(sigCorrPosCells)
    plot(smooth(fineHist(sigCorrPosCells(i),:),11)/max(smooth(fineHist(sigCorrPosCells(i),:),11)) + i/2)
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

%% Find cell types
%now plot things!
[indCellType] = functionCellStringFind(masterHeader,'CellType');


%determine if there are units in each category
% findPVs = find(bigMaster(:,indCellType) == 1);
% 
% findMSNs = find(bigMaster(:,indCellType) == 0);

findCHATs = find(bigMaster(:,indCellType) == 2);


findPVs = find(bigMaster(:,5) < 0.0004 & bigMaster(:,6) > 1.1);
findMSNs = find(bigMaster(:,5) > 0.0005 & bigMaster(:,6) > 1.1); 

%% Pull tuning widths


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
tarDB = 3;
smoothWin = 3;
counter = 1;
for i = 1:length(bigMaster)
    if ismember(i,[1,101,201,301,401,501,601])
        hFig = figure;
        counter = 1;
    end
    subplot(10,10,counter)
%     plotData = squeeze(binValBigStore(:,5,i));
    plotData = binValBigStore(:,:,i);
    plotData = plotData(:,tarDB);
    plotData = smooth(plotData,smoothWin);
    
    if ismember(i,findPVs)
        plot(plotData,'r','LineWidth',2)
    elseif ismember(i,findMSNs)
        plot(plotData,'k','LineWidth',2)
    else
        plot(plotData,'Color',[0.7 0.7 0.7],'LineWidth',2)
    end
    hold on
    plot([1 16],[0 0],'k')
    plot([1 16],[bigWidthCutVal(tarDB,i) bigWidthCutVal(tarDB,i)],'c','LineWidth',2)
    plot([bigWidthMaxPosStore(tarDB,i) - bigWidthSelWidth(tarDB,1,i) bigWidthMaxPosStore(tarDB,i) - bigWidthSelWidth(tarDB,1,i)],[0 max(plotData)],'g','LineWidth',2)
    plot([bigWidthMaxPosStore(tarDB,i) + bigWidthSelWidth(tarDB,2,i) bigWidthMaxPosStore(tarDB,i) + bigWidthSelWidth(tarDB,2,i)],[0 max(plotData)],'g','LineWidth',2)
    plot([bigWidthMaxPosStore(tarDB,i) bigWidthMaxPosStore(tarDB,i)],[0 max(plotData)],'m','LineWidth',2)
    xlim([1 16])
    counter = counter + 1;
end


signStore = sign(binValBigStore);
signSigStore = sigValBigStore.*signStore;
signSigStore(signSigStore <= 0) = NaN;
sigCut = 0.01;
for i = 1:size(signSigStore,3)
    findSigs(i) = length(find(signSigStore(:,:,i) <= sigCut));
end

tarCells = [1:length(bigMaster)];
tarPVs = findPVs;
tarMSNs = findMSNs;
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

%select only for units that have some kind of width at some DB

selPVs = tarPVs;
selPVs(sum(sigWidthPV) == 0) = [];
selMSNs = tarMSNs;
selMSNs(sum(sigWidthMSN) == 0) = [];

selWidthPV = sigWidthPV;
selWidthPV(:,sum(sigWidthPV) == 0) = [];
selWidthMSN = sigWidthMSN;
selWidthMSN(:,sum(sigWidthMSN) == 0) = [];


%hmmm not seeing very significant differences in width here...at amplitude
%= 1 we see something, but thats pretty much it. 


%% Now lets start processing the STA data. 

%find significant dmr stas
zCorr = (bigStoreSTACorr - mean(bigStoreSTACorrBase'))./std(bigStoreSTACorrBase');
for i= 1:length(bigStoreSTACorr)
    corrPrctile(i,:) = prctile(bigStoreSTACorrBase(i,:),[0.5 99.5]);
end

zlim = 3;
sigUnitZ = find(zCorr > zlim);
sigUnitPrct = find(bigStoreSTACorr - corrPrctile(:,2)' > 0);
%these produce more or less overlapping distributions. Can probably use
%iether one. Lets be a bit more lax and use the prct measure. 


findPVs = find(bigMaster(:,5) < 0.0004 & bigMaster(:,6) > 1.1);
findMSNs = find(bigMaster(:,5) > 0.0005 & bigMaster(:,6) > 1.1); 
%first of all, we need to invalidate STAs with too few spikes
minSpikes = 200;
fewSpikes = find(dmrSpikeNum < minSpikes);

dmrMSN = findMSNs;
dmrMSN(ismember(dmrMSN,fewSpikes)) = [];
dmrMSN = intersect(dmrMSN,sigUnitPrct);
dmrPV = findPVs;
dmrPV(ismember(dmrPV,fewSpikes)) = [];
dmrPV = intersect(dmrPV,sigUnitPrct);

%plot out spiking rate for DMR
histVectSpikeNum = [0:5:40];
figure
subplot(2,1,1)
hist(dmrSpikeNum(dmrPV)/(10*60),histVectSpikeNum)
subplot(2,1,2)
hist(dmrSpikeNum(dmrMSN)/(10*60),histVectSpikeNum)

%calculate phase locking index
% PLI = (max(bigSTAstore') - min(bigSTAstore'))./(dmrSpikeNum/(10*60)*sqrt(8));
PLI = (max(bigSTAstore') - min(bigSTAstore'))./(dmrSpikeNum*38.8520);

histVectSpikeNum = [0:0.05:1];
figure
subplot(2,1,1)
hist(PLI(dmrPV),histVectSpikeNum)
subplot(2,1,2)
hist(PLI(dmrMSN),histVectSpikeNum)

%first, reshape the STAs!
for i = 1:length(bigDBStore)
    tmpStore = reshape(bigSTASigstore(i,:),length(faxis),[]);
    newSTA(:,:,i) = tmpStore;
    newTMP = tmpStore;
    newTMP(tmpStore >=0) = 0;
    negSTA(:,:,i) = newTMP;
    tmpStore(tmpStore <= 0) = 0;
    posSTA(:,:,i) = tmpStore;
end

%to determine peak and width, lets compress along time axis. 
smoothWind = 1;
for i = 1:length(bigDBStore)
    compTimingPos(:,i) = sum(posSTA(:,:,i));
    compTimingNeg(:,i) =  sum(negSTA(:,:,i));
    compWidthPos(:,i) = smooth(sum(posSTA(:,:,i)'),smoothWind);
    compWidthNeg(:,i) = smooth(sum(negSTA(:,:,i)'),smoothWind);
    
end
% 
% figure
% hold on
% for i = 1:length(bigDBStore)
%     plot(compWidthPos(:,i))
%     plot(smooth(compWidthPos(:,i),5),'r')
%     
% end
% 
% %peak to peak distance for old datasets seems to be 5. Lets keep a 5
% %smoothing for these. 

%now lets go through the shits. Pull widths and timings. 
dmrWidthPos = [];
dmrBFPos = [];
dmrWidthNeg = [];
dmrBFNeg = [];
widthPer = 0.5;
for i = 1:length(bigDBStore)
    %pull peak modulation timing. 
    [maxVal dmrTimePos(i)] = max(compTimingPos(:,i));
    [minVal dmrTimeNeg(i)] = min(compTimingPos(:,i));
    %pull widths
    [dmrWidthPos(i,1:4),dmrBFPos(i),maxVal,cutVal] = functionHeightBasedTuningWidth(compWidthPos(:,i),widthPer,1);
    [dmrWidthNeg(i,1:4),dmrBFNeg(i),maxVal,cutVal] = functionHeightBasedTuningWidth(-1*compWidthNeg(:,i),widthPer,1);
end

dmrWidthPos (:,[2,4]) = [];
dmrWidthNeg (:,[2,4]) = [];

%first, lets just make plots of tuning curves of pure tones vs DMR.
%start with PVs
axisDef = ceil(sqrt(length(dmrPV)));
hFig = figure;
for i = 1:length(dmrPV)
    %pull tuning curve average from top three amplitudes
    toneWidth = binValBigStore(:,:,dmrPV(i));
    toneWidth = mean(toneWidth(:,end-2:end)');
    %make subplot
    subplot(axisDef,axisDef,i)
    hold on
    plot(toneFreqs,toneWidth/max(toneWidth),'k')
    plot(faxis,compWidthPos(:,dmrPV(i))/max(compWidthPos(:,dmrPV(i))),'r')
    plot(faxis,-1*compWidthNeg(:,dmrPV(i))/min(compWidthNeg(:,dmrPV(i))),'b')
    xlim([faxis(1) faxis(end)])
    set(gca, 'XScale', 'log')
end

%Then MSNs
axisDef = ceil(sqrt(length(dmrMSN)));
hFig = figure;
for i = 1:length(dmrMSN)
    %pull tuning curve average from top three amplitudes
    toneWidth = binValBigStore(:,:,dmrMSN(i));
    toneWidth = mean(toneWidth(:,end-2:end)');
    %make subplot
    subplot(axisDef,axisDef,i)
    hold on
    plot(toneFreqs,toneWidth/max(toneWidth),'k')
    plot(faxis,compWidthPos(:,dmrMSN(i))/max(compWidthPos(:,dmrMSN(i))),'r')
    plot(faxis,-1*compWidthNeg(:,dmrMSN(i))/min(compWidthNeg(:,dmrMSN(i))),'b')
    xlim([faxis(1) faxis(end)])
    set(gca, 'XScale', 'log')
end

%now lets plot out temporal information. 
smoothWind = 6;
timeStep = 9.9800e-04;
%FSI first
axisDef = ceil(sqrt(length(dmrPV)));
hFig = figure;
for i = 1:length(dmrPV)
    %make subplot
    subplot(axisDef,axisDef,i)
    hold on
    plot(fineHistVect,smooth(fineHist(dmrPV(i),:),smoothWind)/max(smooth(fineHist(dmrPV(i),:),smoothWind)),'k')
    plot([timeStep:timeStep:timeStep*100],compTimingPos(end:-1:1,dmrPV(i))/max(compTimingPos(:,dmrPV(i))),'r')
    plot([timeStep:timeStep:timeStep*100],-1*compTimingNeg(end:-1:1,dmrPV(i))/min(compTimingNeg(:,dmrPV(i))),'b')
    xlim([0 0.1])
%     set(gca, 'XScale', 'log')
end

%Then MSNs
axisDef = ceil(sqrt(length(dmrMSN)));
hFig = figure;
for i = 1:length(dmrMSN)
    %make subplot
    subplot(axisDef,axisDef,i)
    hold on
    plot(fineHistVect,smooth(fineHist(dmrMSN(i),:),smoothWind)/max(smooth(fineHist(dmrMSN(i),:),smoothWind)),'k')
    plot([timeStep:timeStep:timeStep*100],compTimingPos(end:-1:1,dmrMSN(i))/max(compTimingPos(:,dmrMSN(i))),'r')
    plot([timeStep:timeStep:timeStep*100],-1*compTimingNeg(end:-1:1,dmrMSN(i))/min(compTimingNeg(:,dmrMSN(i))),'b')
    xlim([0 0.1])
%     set(gca, 'XScale', 'log')
end


%now lets pull out the RTFs
stepSize = timeStep;
timeSteps = 100;
rtf = [];
for i = 1:length(bigDBStore)
    [tmf, xmf, rtf] = sta2rtf(reshape(bigSTAstore(i,:),length(faxis),[]), [stepSize:stepSize:stepSize*100], faxis, 40, 4, 'n');
    rtfStore(:,:,i) = rtf;
    temp1 = flip(rtf(:,1:16),2);
    temp2 = rtf(:,16:end);
    flipRTF(:,:,i) = temp1 + temp2; %second dimension should be temporal mod. First dimension frequency mod. 
    rtfModTemp(:,i) = sum(flipRTF(:,:,i));
    rtfModSpect(:,i) = sum(flipRTF(:,:,i)');
end

%plot out STA, STA SIG, and RTF
%FSI first
% axisDef = ceil(sqrt(length(dmrPV)));

for i = 1:length(dmrPV)
    hFig = figure;
    set(hFig, 'Position', [10 10 500 1000])
    %make subplot
    subplot(3,1,1)
    imagesc(reshape(bigSTAstore(dmrPV(i),:),length(faxis),[]))
    colormap('parula')
    title('STA')
    subplot(3,1,2)
    imagesc(newSTA(:,:,dmrPV(i)))
    colormap('parula')
    title('STA Sig')
    subplot(3,1,3)
    imagesc(flipRTF(:,:,dmrPV(i)))
    colormap('parula')
end

% 
% for i = 1:length(dmrMSN)
%     hFig = figure;
%     set(hFig, 'Position', [10 10 500 1000])
%     %make subplot
%     subplot(3,1,1)
%     imagesc(reshape(bigSTAstore(dmrMSN(i),:),length(faxis),[]))
%     colormap('parula')
%     title('STA')
%     subplot(3,1,2)
%     imagesc(newSTA(:,:,dmrMSN(i)))
%     colormap('parula')
%     title('STA Sig')
%     subplot(3,1,3)
%     imagesc(flipRTF(:,:,dmrMSN(i)))
%     colormap('parula')
% end


axisDef = ceil(sqrt(length(bigDBStore)));
figure
for i = 1:length(bigDBStore)
    subplot(axisDef,axisDef,i)
    if ismember(i,dmrPV)
        plot(rtfModTemp(:,i),'r')
    elseif ismember(i,dmrMSN)
        plot(rtfModTemp(:,i),'k')
    else
        plot(rtfModTemp(:,i),'Color',[0.7 0.7 0.7])
    end
%     plot(rtfModTemp(:,i))
end

axisDef = ceil(sqrt(length(bigDBStore)));
figure
for i = 1:length(bigDBStore)
    subplot(axisDef,axisDef,i)
    if ismember(i,dmrPV)
        plot(rtfModSpect(:,i),'r')
    elseif ismember(i,dmrMSN)
        plot(rtfModSpect(:,i),'k')
    else
        plot(rtfModSpect(:,i),'Color',[0.7 0.7 0.7])
    end
end

figure
axisDef = ceil(sqrt(length(dmrPV)));
for i = 1:length(dmrPV)
    subplot(axisDef,axisDef,i)
    plot(rtfModSpect(:,dmrPV(i)),'r')
end

figure
axisDef = ceil(sqrt(length(dmrMSN)));
for i = 1:length(dmrMSN)
    subplot(axisDef,axisDef,i)
    plot(rtfModSpect(:,dmrMSN(i)),'k')
end


figure
axisDef = ceil(sqrt(length(dmrPV)));
for i = 1:length(dmrPV)
    subplot(axisDef,axisDef,i)
    plot(rtfModTemp(:,dmrPV(i)),'r')
end

figure
axisDef = ceil(sqrt(length(dmrMSN)));
for i = 1:length(dmrMSN)
    subplot(axisDef,axisDef,i)
    plot(rtfModTemp(:,dmrMSN(i)),'k')
end





