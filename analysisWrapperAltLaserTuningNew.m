



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
    
%     
%     %now lets try to get corrcoef from binDiff values (tempBinStore)
%     shuffBinDiffCorrVals = [];
%     binDiffCorrVals = [];
%     binDiffCorrValSig = [];
%     for j = 1:numUnits
%         val1 = reshape(squeeze(tempBinStore(:,:,j)),1,[]);
%         for k = 1:numUnits
%             val2 = reshape(squeeze(tempBinStore(:,:,k)),1,[]);
%             tempCorrVal = corrcoef(val1,val2);
%             binDiffCorrVals(j,k) = tempCorrVal(2);
%             %now do shuffled versions
%             numShuff = 1000;
%             for l = 1:numShuff
%                 shuffVal = randperm(length(val2));
%                 shuffVal = val2(shuffVal);
%                 tempCorrVal = corrcoef(val1,shuffVal);
%                 shuffBinDiffCorrVals(l) = tempCorrVal(2);
%             end
%             %calculate percentile
%             prctileVals = prctile(shuffBinDiffCorrVals,[0.5 99.5]);
%             if binDiffCorrVals(j,k) < prctileVals(1) || binDiffCorrVals(j,k) > prctileVals(2)
%                 binDiffCorrValSig(j,k) = 1;
%             else
%                 binDiffCorrValSig(j,k) = 0;
%             end
%             shuffBinDiffCorrVals = [];
%             %also calculate distances between units
%             distVal(j,k) = sqrt(((floor(masterData(j,1)) - floor(masterData(k,1)))*25)^2 + ((masterData(j,2) - masterData(k,2))*250)^2);
%         end
%     end
%     findPVs = find(masterData(:,5) < 0.0004 & masterData(:,6) > 1.1);
%     findMSNs = find(masterData(:,5) > 0.0005 & masterData(:,6) > 1.1); 
%     findCHATs = find(masterData(:,6) < 1.1);
%     corrCoefData.CorrCoefs = binDiffCorrVals;
%     corrCoefData.CorrCoefSig = binDiffCorrValSig;
%     corrCoefData.CellType = NaN(numUnits,1);
%     corrCoefData.CellType(findPVs) = 1;
%     corrCoefData.CellType(findMSNs) = 0;
%     corrCoefData.CellType(findCHATs) = 4;
%     corrCoefData.Distance = distVal;
%     bigCorrStore{i} = corrCoefData;
%     corrCoefData = [];
%     distVal = [];
    
    bigMasterInd = bigMasterInd + numUnits;
end

findMSNs = find(bigMaster(:,7) == 0);
findPVs = find(bigMaster(:,7) == 1);


%% now lets march through the correlation coefficient data and see what we pull out. First things first, I want to remove all the excess values, since I dont want double counting. 
%lets define interactions! There are three cell types, and therefore six
%interactions: FSI-MSN, FSI-ChAT, FSI-FSI, MSN-ChAT, MSN-MSN, ChAT-ChAT. 
% 
% %lets number these. I want the most common first, going up into the less
% %common. so:
% 
% % MSN-MSN: 0
% % MSN-FSI: 1
% % FSI-FSI: 2
% % MSN-ChAT: 4
% % FSI-ChAT: 5
% % ChAT-ChAT: 8
% 
% 
% corrCount = 1;
% cellCount = 1;
% for i = 1:length(bigCorrStore)
%     dataset = bigCorrStore{i};
%     %now lets try and linearize this dataset. 
%     numUnits = length(dataset.CellType);
%     for j = 1:numUnits
%         corrValStore(corrCount:corrCount + numUnits - j - 1) = dataset.CorrCoefs(j,j+1:end);
%         distValStore(corrCount:corrCount + numUnits - j - 1) = dataset.Distance(j,j+1:end);
%         interTypeStore(corrCount:corrCount + numUnits - j - 1) = dataset.CellType(j)+ dataset.CellType(j+1:end);
%         corrValSigStore(corrCount:corrCount + numUnits - j - 1) = dataset.CorrCoefSig(j,j+1:end);
%         cellTrackerStore(corrCount:corrCount + numUnits - j - 1,1) = i;
%         cellTrackerStore(corrCount:corrCount + numUnits - j - 1,2) = j+cellCount-1;
%         cellTrackerStore(corrCount:corrCount + numUnits - j - 1,3) = cellCount + (j+1:numUnits) - 1;
%         cellTrackerStore(corrCount:corrCount + numUnits - j - 1,4) = dataset.CellType(j);
%         cellTrackerStore(corrCount:corrCount + numUnits - j - 1,5) = dataset.CellType(j+1:numUnits);
%         cellTrackerStore(corrCount:corrCount + numUnits - j - 1,6) = bigMaster(j+cellCount-1,2);
%         cellTrackerStore(corrCount:corrCount + numUnits - j - 1,7) = bigMaster(cellCount + (j+1:numUnits) - 1,2);
%         
%         corrCount = corrCount + numUnits - j;
%     end
%     cellCount = cellCount + numUnits;
% end
% 
% %now lets clean up the dataset. Remove all NaNs. 
% corrValStore(isnan(interTypeStore)) = [];
% distValStore(isnan(interTypeStore)) = [];
% cellTrackerStore(isnan(interTypeStore),:) = [];
% corrValSigStore(isnan(interTypeStore)) = [];
% interTypeStore(isnan(interTypeStore)) = [];
% 
% %lets look in general at within neuron corrcoefs. 
% corrHistVect = [-1:0.05:1];
% 
% %just MSNs, close
% selDist = 150;
% selInterType = 0;
% 
% findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
% figure
% corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);
% sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
% corrValmsnmsnSig = hist(corrValStore(sigInt),corrHistVect);
% hold on
% bar(corrHistVect,corrValmsnmsn)
% bar(corrHistVect,corrValmsnmsnSig,'k')
% % hist(corrValStore(findTarCorr),corrHistVect)
% xlim([-1 1])
% % corrValmsnmsn = hist(corrValStore(findTarCorr),corrHistVect);
% 
% %just FSIs, close
% selDist = 150;
% selInterType = 2;
% 
% 
% findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
% figure
% corrValfsifsi = hist(corrValStore(findTarCorr),corrHistVect);
% sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
% corrValfsifsiSig = hist(corrValStore(sigInt),corrHistVect);
% hold on
% bar(corrHistVect,corrValfsifsi)
% bar(corrHistVect,corrValfsifsiSig,'k')
% % hist(corrValStore(findTarCorr),corrHistVect)
% xlim([-1 1])
% 
% 
% %MSN vs FSI, close
% selDist = 150;
% selInterType = 1;
% 
% 
% findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
% figure
% corrValfsimsn = hist(corrValStore(findTarCorr),corrHistVect);
% sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
% corrValfsimsnSig = hist(corrValStore(sigInt),corrHistVect);
% hold on
% bar(corrHistVect,corrValfsimsn)
% bar(corrHistVect,corrValfsimsnSig,'k')
% % hist(corrValStore(findTarCorr),corrHistVect)
% xlim([-1 1])
% 
% %MSN vs ChAT, close
% selDist = 150;
% selInterType = 5;
% 
% 
% findTarCorr = intersect(find(distValStore <= selDist),find(interTypeStore == selInterType));
% figure
% corrValmsnchat = hist(corrValStore(findTarCorr),corrHistVect);
% sigInt = intersect(findTarCorr,find(corrValSigStore == 1));
% corrValmsnchatSig = hist(corrValStore(sigInt),corrHistVect);
% hold on
% bar(corrHistVect,corrValmsnchat)
% bar(corrHistVect,corrValmsnchatSig,'k')
% xlim([-1 1])
%     
% 
% %plot cumdists?
% figure
% hold on
% plot(corrHistVect,cumsum(corrValmsnmsn)/sum(corrValmsnmsn),'k')
% plot(corrHistVect,cumsum(corrValfsifsi)/sum(corrValfsifsi),'c')
% plot(corrHistVect,cumsum(corrValfsimsn)/sum(corrValfsimsn),'r')
% 
% %try plotting by distance for first 150. 
% selInterType = 1;
% for i = 1:9
%     selDist = (i-1)*25;
%     tempTar = intersect(find(distValStore == selDist),find(interTypeStore == selInterType));
%     distCorrHistStore(:,i) = hist(corrValStore(tempTar),corrHistVect);
%     normDistCorrHistStore(:,i) = distCorrHistStore(:,i)/max(distCorrHistStore(:,i));
% end
% 
% figure
% for i = 1:9
%     subplot(9,1,i)
%     plot(corrHistVect,normDistCorrHistStore(:,i));
% end
% 
% %now lets try and dig into the units that I'm pulling out. Lets pull out
% %significant MSN-FSI and split between negative and positive. 
% selInterType = 1;
% selDist = 200;
% findSigCorrNeg = find(corrValSigStore == 1 & interTypeStore == selInterType & corrValStore < 0 & distValStore < selDist);
% findSigCorrPos = find(corrValSigStore == 1 & interTypeStore == selInterType & corrValStore > 0 & distValStore < selDist);
% 
% 
% %now we need to look at cellTrackerStore and extract the MSNs
% for i = 1:length(findSigCorrNeg)
%     testCells = cellTrackerStore(findSigCorrNeg(i),:);
%     if testCells(4) == 0
%         sigCorrNegCells(i) = testCells(2);
%     elseif testCells(5) == 0
%         sigCorrNegCells(i) = testCells(3);
%     end
% end
% 
% for i = 1:length(findSigCorrPos)
%     testCells = cellTrackerStore(findSigCorrPos(i),:);
%     if testCells(4) == 0
%         sigCorrPosCells(i) = testCells(2);
%     elseif testCells(5) == 0
%         sigCorrPosCells(i) = testCells(3);
%     end
% end
% %this process can generate duplicates if a single cell has multiple
% %negative correlations with neighboring FSIs. Eliminate duplicates. 
% sigCorrNegCells = unique(sigCorrNegCells);
% sigCorrPosCells = unique(sigCorrPosCells);
% 
% 
% %now lets plot out tuning curves
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.03 0.03], [0.03 0.03]);
% axisSize = ceil(sqrt(length(sigCorrNegCells)));
% for i = 1:length(sigCorrNegCells)
%     subplot(axisSize,axisSize,i)
%     imagesc(binValBigStore{sigCorrNegCells(i)}')
%     colormap('parula')
% end
% 
% hFig = figure;
% axisSize = ceil(sqrt(length(sigCorrPosCells)));
% for i = 1:length(sigCorrPosCells)
%     subplot(axisSize,axisSize,i)
%     imagesc(binValBigStore{sigCorrPosCells(i)}')
%     colormap('parula')
% end
% 
% %now lets look at the MSN-MSN interactions of the targeted cells. 
% selDist = 200;
% selType = 0;
% negCellmsnCorr = [];
% negCellmsnCorrSig = [];
% counter = 1;
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% for i = 1:length(sigCorrNegCells)
%     findRow1 = find(cellTrackerStore(:,2) == sigCorrNegCells(i));
%     findRow2 = find(cellTrackerStore(:,3) == sigCorrNegCells(i));
%     
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% counter = 1;
% posCellmsnCorr = [];
% posCellmsnCorrSig = [];
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% for i = 1:length(sigCorrPosCells)
%     findRow1 = find(cellTrackerStore(:,2) == sigCorrPosCells(i));
%     findRow2 = find(cellTrackerStore(:,3) == sigCorrPosCells(i));
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% hFig = figure;
% hold on
% negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
% negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,negCellBarPlot)
% bar(corrHistVect,negCellSigBarPlot,'k')
% xlim([-1 1])
% 
% hFig = figure;
% hold on
% posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
% posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,posCellBarPlot)
% bar(corrHistVect,posCellSigBarPlot,'k')
% xlim([-1 1])
% 
% %now lets do the same, but eliminate possibility for double counting
% %interactions, which may exist now. 
% selDist = 200;
% selType = 0;
% negCellmsnCorr = [];
% negCellmsnCorrSig = [];
% counter = 1;
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% tempTrack = cellTrackerStore;
% for i = 1:length(sigCorrNegCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
%     findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% counter = 1;
% posCellmsnCorr = [];
% posCellmsnCorrSig = [];
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% tempTrack = cellTrackerStore;
% for i = 1:length(sigCorrPosCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
%     findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% hFig = figure;
% hold on
% negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
% negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,negCellBarPlot)
% bar(corrHistVect,negCellSigBarPlot,'k')
% xlim([-1 1])
% 
% hFig = figure;
% hold on
% posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
% posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,posCellBarPlot)
% bar(corrHistVect,posCellSigBarPlot,'k')
% xlim([-1 1])
% 
% %now lets try and do it only for cells within the same group. 
% selDist = 200;
% negCellmsnCorr = [];
% negCellmsnCorrSig = [];
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% tempTrack = cellTrackerStore;
% counter = 1;
% for i = 1:length(sigCorrNegCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
%     if findRow1
%         otherVals = tempTrack(findRow1,3);
%         otherVals = find(ismember(otherVals,sigCorrNegCells));
%         findRow1 = findRow1(otherVals);
%     end
%     
%     
%     findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
%     if findRow2
%         otherVals = tempTrack(findRow2,2);
%         otherVals = find(ismember(otherVals,sigCorrNegCells));
%         findRow2 = findRow2(otherVals);
%     end
%     
%     
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% counter = 1;
% posCellmsnCorr = [];
% posCellmsnCorrSig = [];
% tempTrack = cellTrackerStore;
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% for i = 1:length(sigCorrPosCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
%     if findRow1
%         otherVals = tempTrack(findRow1,3);
%         otherVals = find(ismember(otherVals,sigCorrPosCells));
%         findRow1 = findRow1(otherVals);
%     end
%     
%     
%     findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
%     if findRow2
%         otherVals = tempTrack(findRow2,2);
%         otherVals = find(ismember(otherVals,sigCorrPosCells));
%         findRow2 = findRow2(otherVals);
%     end
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% hFig = figure;
% hold on
% negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
% negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,negCellBarPlot)
% bar(corrHistVect,negCellSigBarPlot,'k')
% xlim([-1 1])
% 
% hFig = figure;
% hold on
% posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
% posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,posCellBarPlot)
% bar(corrHistVect,posCellSigBarPlot,'k')
% xlim([-1 1])
% 
% %now lets try and do it only for cells not in the same group. 
% selDist = 200;
% negCellmsnCorr = [];
% negCellmsnCorrSig = [];
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% tempTrack = cellTrackerStore;
% counter = 1;
% for i = 1:length(sigCorrNegCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
%     if findRow1
%         otherVals = tempTrack(findRow1,3);
%         otherVals = find(ismember(otherVals,sigCorrPosCells));
%         findRow1 = findRow1(otherVals);
%     end
%     
%     
%     findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
%     if findRow2
%         otherVals = tempTrack(findRow2,2);
%         otherVals = find(ismember(otherVals,sigCorrPosCells));
%         findRow2 = findRow2(otherVals);
%     end
%     
%     
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% counter = 1;
% posCellmsnCorr = [];
% posCellmsnCorrSig = [];
% tempTrack = cellTrackerStore;
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% for i = 1:length(sigCorrPosCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
%     if findRow1
%         otherVals = tempTrack(findRow1,3);
%         otherVals = find(ismember(otherVals,sigCorrNegCells));
%         findRow1 = findRow1(otherVals);
%     end
%     
%     
%     findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
%     if findRow2
%         otherVals = tempTrack(findRow2,2);
%         otherVals = find(ismember(otherVals,sigCorrNegCells));
%         findRow2 = findRow2(otherVals);
%     end
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% hFig = figure;
% hold on
% negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
% negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,negCellBarPlot)
% bar(corrHistVect,negCellSigBarPlot,'k')
% xlim([-1 1])
% 
% hFig = figure;
% hold on
% posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
% posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,posCellBarPlot)
% bar(corrHistVect,posCellSigBarPlot,'k')
% xlim([-1 1])
% 
% %now lets try and do it only for cells not in either group. 
% selDist = 200;
% negCellmsnCorr = [];
% negCellmsnCorrSig = [];
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% tempTrack = cellTrackerStore;
% counter = 1;
% for i = 1:length(sigCorrNegCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrNegCells(i));
%     if findRow1
%         otherVals = tempTrack(findRow1,3);
%         otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
%         findRow1 = findRow1(otherVals);
%     end
%     
%     
%     findRow2 = find(tempTrack(:,3) == sigCorrNegCells(i));
%     if findRow2
%         otherVals = tempTrack(findRow2,2);
%         otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
%         findRow2 = findRow2(otherVals);
%     end
%     
%     
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     negCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     negCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% counter = 1;
% posCellmsnCorr = [];
% posCellmsnCorrSig = [];
% tempTrack = cellTrackerStore;
% findRow3 = find(distValStore < selDist)';
% findRow4 = find(interTypeStore == selType)';
% findRow3 = intersect(findRow3,findRow4);
% for i = 1:length(sigCorrPosCells)
%     findRow1 = find(tempTrack(:,2) == sigCorrPosCells(i));
%     if findRow1
%         otherVals = tempTrack(findRow1,3);
%         otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
%         findRow1 = findRow1(otherVals);
%     end
%     
%     
%     findRow2 = find(tempTrack(:,3) == sigCorrPosCells(i));
%     if findRow2
%         otherVals = tempTrack(findRow2,2);
%         otherVals = find(~ismember(otherVals,sigCorrPosCells) & ~ismember(otherVals,sigCorrNegCells));
%         findRow2 = findRow2(otherVals);
%     end
% %     findRow3 = find(distValStore < selDist)';
%     allFinds = sort([findRow1;findRow2]);
%     allFinds = intersect(allFinds,findRow3);
%     tempTrack(allFinds,:) = 0;
%     posCellmsnCorr(counter:counter + length(allFinds) - 1) = corrValStore(allFinds);
%     posCellmsnCorrSig(counter:counter + length(allFinds) - 1) = corrValSigStore(allFinds);
%     counter = counter + length(allFinds);
% end
% 
% hFig = figure;
% hold on
% negCellBarPlot = hist(negCellmsnCorr,corrHistVect);
% negCellSigBarPlot = hist(negCellmsnCorr(negCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,negCellBarPlot)
% bar(corrHistVect,negCellSigBarPlot,'k')
% xlim([-1 1])
% 
% hFig = figure;
% hold on
% posCellBarPlot = hist(posCellmsnCorr,corrHistVect);
% posCellSigBarPlot = hist(posCellmsnCorr(posCellmsnCorrSig==1),corrHistVect);
% bar(corrHistVect,posCellBarPlot)
% bar(corrHistVect,posCellSigBarPlot,'k')
% xlim([-1 1])
% 
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

%now lets try and plot out things!
% tarDB = 3;
% smoothWin = 3;
% counter = 1;
% for i = 1:length(bigMaster)
%     if ismember(i,[1,101,201,301,401,501,601])
%         hFig = figure;
%         counter = 1;
%     end
%     subplot(10,10,counter)
% %     plotData = squeeze(binValBigStore(:,5,i));
%     plotData = binValBigStore{i};
%     plotData = plotData(:,tarDB);
%     plotData = smooth(plotData,smoothWin);
%     
%     if ismember(i,findPVs)
%         plot(plotData,'r','LineWidth',2)
%     elseif ismember(i,findMSNs)
%         plot(plotData,'k','LineWidth',2)
%     else
%         plot(plotData,'Color',[0.7 0.7 0.7],'LineWidth',2)
%     end
%     hold on
%     plot([1 freqNumStore(i)],[0 0],'k')
%     plot([1 freqNumStore(i)],[bigWidthCutVal(tarDB,i) bigWidthCutVal(tarDB,i)],'c','LineWidth',2)
%     plot([bigWidthMaxPosStore(tarDB,i) - bigWidthSelWidth(tarDB,1,i) bigWidthMaxPosStore(tarDB,i) - bigWidthSelWidth(tarDB,1,i)],[0 max(plotData)],'g','LineWidth',2)
%     plot([bigWidthMaxPosStore(tarDB,i) + bigWidthSelWidth(tarDB,2,i) bigWidthMaxPosStore(tarDB,i) + bigWidthSelWidth(tarDB,2,i)],[0 max(plotData)],'g','LineWidth',2)
%     plot([bigWidthMaxPosStore(tarDB,i) bigWidthMaxPosStore(tarDB,i)],[0 max(plotData)],'m','LineWidth',2)
%     xlim([1 freqNumStore(i)])
%     counter = counter + 1;
% end
% 
% 
% %now lets try and plot out things!
% tarDB = 1;
% smoothWin = 3;
% counter = 1;
% for i = 1:length(bigMaster)
%     if ismember(i,[1,101,201,301,401,501,601])
%         hFig = figure;
%         counter = 1;
%     end
%     subplot(10,10,counter)
% %     plotData = squeeze(binValBigStore(:,5,i));
%     plotData = binValBigStoreLaser{i};
%     plotData = plotData(:,tarDB);
%     plotData = smooth(plotData,smoothWin);
%     
%     if ismember(i,findPVs)
%         plot(plotData,'r','LineWidth',2)
%     elseif ismember(i,findMSNs)
%         plot(plotData,'k','LineWidth',2)
%     else
%         plot(plotData,'Color',[0.7 0.7 0.7],'LineWidth',2)
%     end
%     hold on
%     plot([1 freqNumStore(i)],[0 0],'k')
%     plot([1 freqNumStore(i)],[bigWidthCutVal(tarDB,i) bigWidthCutVal(tarDB,i)],'c','LineWidth',2)
%     plot([bigWidthMaxPosStoreLaser(tarDB,i) - bigWidthSelWidthLaser(tarDB,1,i) bigWidthMaxPosStoreLaser(tarDB,i) - bigWidthSelWidthLaser(tarDB,1,i)],[0 max(plotData)],'g','LineWidth',2)
%     plot([bigWidthMaxPosStoreLaser(tarDB,i) + bigWidthSelWidthLaser(tarDB,2,i) bigWidthMaxPosStoreLaser(tarDB,i) + bigWidthSelWidthLaser(tarDB,2,i)],[0 max(plotData)],'g','LineWidth',2)
%     plot([bigWidthMaxPosStoreLaser(tarDB,i) bigWidthMaxPosStoreLaser(tarDB,i)],[0 max(plotData)],'m','LineWidth',2)
%     xlim([1 freqNumStore(i)])
%     counter = counter + 1;
% end


%% Lets look at width values

%now lets determine whether significant values exist at each DB level.
signStore= [];
signSig = [];
signSigStore = [];
sigCut = 0.05;
for i = 1:length(binValBigStore)
    signStore = sign(binValBigStore{i});
    signSig = signStore.* sigValBigStore{i};
    signSig(signSig <= 0) = NaN;
    signSigStore(:,i) = min(signSig);
    signSigBinaryStore(:,i) = min(signSig);
    signSigBinaryStore(min(signSig) < sigCut,i) = 10;
    signSigBinaryStore(signSigBinaryStore(:,i) < 10,i) = 0;
    signSigBinaryStore(signSigBinaryStore(:,i) == 10,i) = 1;
end
signSigBinaryStore(isnan(signSigBinaryStore)) = 0;

%ok, now we have significant values, max points, and the detected width in
%each direction. Lets try and pull things out now?

sigWidthStore = zeros(size(signSigBinaryStore));
sideWarnLow = zeros(size(signSigBinaryStore));
sideWarnHi = zeros(size(signSigBinaryStore));
for i = 1:length(binValBigStore)
    for j = 1:size(signSigBinaryStore,1)
        %first, check if there is a significant value
        if signSigBinaryStore(j,i) == 1
            %now check for width values. If NaN on either side, determine
            %edge.
            %start with lower edge. 
            if isnan(bigWidthSelWidth(j,1,i))
                lowWidth = bigWidthMaxPosStore(j,i);
                sideWarnLow(j,i) = 1;
            else
                lowWidth = bigWidthSelWidth(j,1,i);
            end
            if isnan(bigWidthSelWidth(j,2,i))
                hiWidth = 16 - bigWidthMaxPosStore(j,i);
                sideWarnHi(j,i) = 1;
            else
                hiWidth = bigWidthSelWidth(j,2,i);
            end
            sigWidthStore(j,i) = hiWidth + lowWidth;
        end
    end
end

%now separate to MSNs vs FSIs
sigWidthPV = sigWidthStore(:,findPVs);
sigWidthMSN = sigWidthStore(:,findMSNs);

%now lets just do some ugly stuff and look only at 70dB point. 
widthValsPV = sigWidthPV(end,:);
widthValsMSN = sigWidthMSN(end,:);

widthHistVect = [0:0.5:16];
histValPV = hist(widthValsPV,widthHistVect);
histValMSN = hist(widthValsMSN,widthHistVect);

%do the same for laser

%now lets determine whether significant values exist at each DB level.
signStore= [];
signSig = [];
signSigStore = [];
sigCut = 0.05;
for i = 1:length(binValBigStoreLaser)
    signStore = sign(binValBigStoreLaser{i});
    signSig = signStore.* sigValBigStoreLaser{i};
    signSig(signSig <= 0) = NaN;
    signSigStore(:,i) = min(signSig);
    signSigBinaryStoreLaser(:,i) = min(signSig);
    signSigBinaryStoreLaser(min(signSig) < sigCut,i) = 10;
    signSigBinaryStoreLaser(signSigBinaryStoreLaser(:,i) < 10,i) = 0;
    signSigBinaryStoreLaser(signSigBinaryStoreLaser(:,i) == 10,i) = 1;
end
signSigBinaryStoreLaser(isnan(signSigBinaryStoreLaser)) = 0;

%ok, now we have significant values, max points, and the detected width in
%each direction. Lets try and pull things out now?

sigWidthStoreLaser = zeros(size(signSigBinaryStoreLaser));
sideWarnLow = zeros(size(signSigBinaryStoreLaser));
sideWarnHi = zeros(size(signSigBinaryStoreLaser));
for i = 1:length(binValBigStoreLaser)
    for j = 1:size(signSigBinaryStoreLaser,1)
        %first, check if there is a significant value
        if signSigBinaryStoreLaser(j,i) == 1
            %now check for width values. If NaN on either side, determine
            %edge.
            %start with lower edge. 
            if isnan(bigWidthSelWidthLaser(j,1,i))
                lowWidth = bigWidthMaxPosStoreLaser(j,i);
                sideWarnLow(j,i) = 1;
            else
                lowWidth = bigWidthSelWidthLaser(j,1,i);
            end
            if isnan(bigWidthSelWidthLaser(j,2,i))
                hiWidth = 16 - bigWidthMaxPosStoreLaser(j,i);
                sideWarnHi(j,i) = 1;
            else
                hiWidth = bigWidthSelWidthLaser(j,2,i);
            end
            sigWidthStoreLaser(j,i) = hiWidth + lowWidth;
        end
    end
end

%now separate to MSNs vs FSIs
sigWidthPVLaser = sigWidthStoreLaser(:,findPVs);
sigWidthMSNLaser = sigWidthStoreLaser(:,findMSNs);

%now lets just do some ugly stuff and look only at 70dB point. 
widthValsPVLaser = sigWidthPVLaser(end,:);
widthValsMSNLaser = sigWidthMSNLaser(end,:);

widthHistVect = [0:0.5:16];
histValPVLaser = hist(widthValsPVLaser,widthHistVect);
histValMSNLaser = hist(widthValsMSNLaser,widthHistVect);


% signrank(sigWidthMSN(1,:),sigWidthMSNLaser(1,:))
% signrank(sigWidthMSN(2,:),sigWidthMSNLaser(2,:))
% signrank(sigWidthMSN(3,:),sigWidthMSNLaser(3,:))

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

% tarCells = find(numSig > 5);
% tarPVs = intersect(tarCells,findPVs);
% tarMSNs = intersect(tarCells,findMSNs);
% 
% signrank(sigWidthStore(1,tarMSNs),sigWidthStoreLaser(1,tarMSNs))
% signrank(sigWidthStore(2,tarMSNs),sigWidthStoreLaser(2,tarMSNs))
% signrank(sigWidthStore(3,tarMSNs),sigWidthStoreLaser(3,tarMSNs))
% 
% signrank(sigWidthStore(1,tarPVs),sigWidthStoreLaser(1,tarPVs))
% signrank(sigWidthStore(2,tarPVs),sigWidthStoreLaser(2,tarPVs))
% signrank(sigWidthStore(3,tarPVs),sigWidthStoreLaser(3,tarPVs))

%NOT SIGNIFICANT!


%now lets go through and see which units had significant changes to size of
%response
for i = 1:length(binValBigStoreLaser)
    sigAmpStore(i,1) = signrank(reshape(binValBigStoreLaser{i},1,[]),reshape(binValBigStore{i},1,[]));
    testData = binValBigStoreLaser{i} - binValBigStore{i};
    changeDir(i) = sign(mean(mean(testData)));
    sigAmpStore(i,2) = signrank(reshape(testData,1,[]));
    
end

tarCells = find(sigAmpStore(:,1) < 0.05);
tarPVs = intersect(tarCells,findPVs);
tarMSNs = intersect(tarCells,findMSNs);

signrank(sigWidthStore(1,tarMSNs),sigWidthStoreLaser(1,tarMSNs))
signrank(sigWidthStore(2,tarMSNs),sigWidthStoreLaser(2,tarMSNs))
signrank(sigWidthStore(3,tarMSNs),sigWidthStoreLaser(3,tarMSNs))

%plot out selected time histograms. 
for i = 1:length(tarMSNs)
figure
subplot(2,1,1)
hold on
plot(smooth(fineHist(tarMSNs(i),:),11),'k')
plot(smooth(fineHistLaser(tarMSNs(i),:),11),'g')
subplot(2,1,2)
plot(smooth(fineHistLaser(tarMSNs(i),:) - fineHist(tarMSNs(i),:),11))
end

%plot out all MSNs subtraction
figure
for i = 1:length(findMSNs)
    subplot(14,14,i)
    plot(smooth(fineHistLaser(findMSNs(i),:) - fineHist(findMSNs(i),:),21))
end



figure
hold on
tarDB = 3;
for i = 1:length(tarMSNs)
    plot([1,2],[sigWidthStore(tarDB,tarMSNs(i)) sigWidthStoreLaser(tarDB,tarMSNs(i))])
end


tarDB = 2;
for i = 1:3
    tarDB = i;
    finder = find(sigWidthStore(tarDB,tarMSNs) == 0);
    find2 = find(sigWidthStoreLaser(tarDB,tarMSNs) == 0);
    finder = unique([finder,find2]);
    tester = sigWidthStore(:,tarMSNs);
    test2 = sigWidthStoreLaser(:,tarMSNs);
    tester(:,finder) = [];
    test2(:,finder) = [];
    signrank(tester(tarDB,:),test2(tarDB,:))
    hFig = figure;
    plot(tester(tarDB,:),test2(tarDB,:),'r.')
    hold on
    plot([0 max(tester(tarDB,:))],[0 max(tester(tarDB,:))],'k')
    spikeGraphName = strcat('ScatterWidthFor',num2str(i));
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end



signrank(sigWidthStore(1,tarPVs),sigWidthStoreLaser(1,tarPVs))
signrank(sigWidthStore(2,tarPVs),sigWidthStoreLaser(2,tarPVs))
signrank(sigWidthStore(3,tarPVs),sigWidthStoreLaser(3,tarPVs))

tarDB = 1;
figure
hold on
for i = 1:tarMSNs
plot([1 2],[0 sigWidthStoreLaser(tarDB,tarMSNs(i))-sigWidthStore(tarDB,tarMSNs(i))])
end
plot([1 2],[0 0],'k','LineWidth',2)


figure
hold on
for i = 1:length(findMSNs)
plot([1,2],[sigWidthMSN(1,i) sigWidthMSNLaser(1,i)])
end
figure
hold on
for i = 1:length(findMSNs)
plot([1,2],[sigWidthMSN(2,i) sigWidthMSNLaser(2,i)])
end
figure
hold on
for i = 1:length(findMSNs)
plot([1,2],[sigWidthMSN(3,i) sigWidthMSNLaser(3,i)])
end

tarVol = 1;
tester = find(sigWidthMSN(tarVol,:) < 1);
tester2 = find(sigWidthMSNLaser(tarVol,:) < 1);
tester = unique([tester,tester2]);
testWidth = sigWidthMSN(tarVol,:);
testWidth(tester) = [];
testWidthLaser = sigWidthMSNLaser(tarVol,:);
testWidthLaser(tester) = [];
mean(testWidthLaser) - mean(testWidth)
signrank(testWidth,testWidthLaser)

figure
plot(sigWidthMSN(1,:),sigWidthMSNLaser(1,:),'k.')
hold on
plot(sigWidthMSN(2,:),sigWidthMSNLaser(2,:),'c.')
plot(sigWidthMSN(3,:),sigWidthMSNLaser(3,:),'m.')

tester = find(sigWidthMSN(1,:) == 0);
test2 = find(sigWidthMSNLaser(1,:) > 0);
growCells = intersect(tester,test2);



%% Now lets try looking at aligning tuning curves


%BF readings from master data are very suspect. Instead, lets go through
%binValBigStore

bigStoreMSN = NaN(length(tarMSNs),41);
bigStoreLaserMSN = NaN(length(tarMSNs),41);

tarDB = 2;
for i = 1:length(tarMSNs)
    %pull curves
    tmp1 = binValBigStore{tarMSNs(i)}(:,tarDB);
    tmp2 = binValBigStoreLaser{tarMSNs(i)}(:,tarDB);
    %find relative maxes. 
    [maxVal1 maxInd1] = max(tmp1);
    [maxVal2 maxInd2] = max(tmp2);
    [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
    bfSave(i,:) = [maxInd1,maxInd2];
    %divide by maximum value
    if maxVal1 > 0
        bigStoreMSN(i,21-maxInd1+1:21-maxInd1+length(tmp1)) = tmp1/maxVal1;
        bigStoreLaserMSN(i,21-maxInd1+1:21-maxInd1+length(tmp1)) = tmp2/maxVal1;
    end
    
end

figure
plot(nanmean(bigStoreMSN))
hold on
plot(nanmean(bigStoreLaserMSN),'r')


%target each one separately


bigStoreMSN = NaN(length(tarMSNs),41);
bigStoreLaserMSN = NaN(length(tarMSNs),41);

% tarDB = 2;
for i = 1:length(tarMSNs)
    %pull curves
    tmp1 = binValBigStore{tarMSNs(i)}(:,tarDB);
    tmp2 = binValBigStoreLaser{tarMSNs(i)}(:,tarDB);
    %find relative maxes. 
    [maxVal1 maxInd1] = max(tmp1);
    [maxVal2 maxInd2] = max(tmp2);
    [maxmaxVal maxmaxInd] = max([maxVal1 maxVal2]);
    
    %divide by maximum value
    if maxVal1 > 0
        bigStoreMSN(i,21-maxInd1+1:21-maxInd1+length(tmp1)) = tmp1/maxVal1;
        
    end
    if maxVal2 > 0
        bigStoreLaserMSN(i,21-maxInd2+1:21-maxInd2+length(tmp1)) = tmp2/maxVal2;
    end
end

figure
plot(nanmean(bigStoreMSN))
hold on
plot(nanmean(bigStoreLaserMSN),'r')




