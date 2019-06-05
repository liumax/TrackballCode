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
        %store number of spikes total
        spikeNum(bigMasterInd + j - 1) = length(s.(desigName{j}).SpikeTimes);
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
        crudeHist(bigMasterInd + j - 1,:) = hist(tempFineRast(:,1),[-0.2:0.005:0.4])/0.005/length(toneFinder);
        nameStore{bigMasterInd + j -1} = targetFiles{i};
        recNumStore(bigMasterInd+j-1) = i;
        unitStore{bigMasterInd + j -1} = desigName{j};
    end
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

%% Lets look at the waveforms a bit more carefully to determine cell types

load('PeakTroughAndHalfWidthComp.mat') %column 1 is from filtered data, column 2 is from unfiltered broadband

[indCellType] = functionCellStringFind(masterHeader,'CellType');


%determine if there are units in each category
% findPVs = find(bigMaster(:,indCellType) == 1);
% 
% findMSNs = find(bigMaster(:,indCellType) == 0);

findCHATs = find(peakTrough(:,2) > 0.0005 & bigMaster(:,6) < 1.1);


findPVs = find(peakTrough(:,2) < 0.00055 & bigMaster(:,6) > 1.1);
findMSNs = find(peakTrough(:,2) > 0.0006 & bigMaster(:,6) > 1.1); 

%eliminate off bad RPVs
badRPV = find(rpvStore > 0.5);
findBads = intersect(findPVs,badRPV);
findPVs(ismember(findPVs,findBads)) = [];
findBads = intersect(findMSNs,badRPV);
findMSNs(ismember(findMSNs,findBads)) = [];

%also eliminate off bad number of spikes. 
minSpikeNum = 250;
badSpikes = find(spikeNum < minSpikeNum);
findBads = intersect(findPVs,badSpikes);
findPVs(ismember(findPVs,findBads)) = [];
findBads = intersect(findMSNs,badSpikes);
findMSNs(ismember(findMSNs,findBads)) = [];

%% Now really crude look: Lets just make population PSTHs.

popPSTHfsi = mean(fineHist(findPVs,:));
popPSTHmsn = mean(fineHist(findMSNs,:));

hFig = figure;
subplot(2,1,1)
hold on
plot([-0.2:binSize:0.4],popPSTHfsi,'r','LineWidth',2)
plot([-0.2:binSize:0.4],popPSTHfsi + std(fineHist(findPVs,:))/sqrt(length(findPVs)),'r','LineWidth',1)
plot([-0.2:binSize:0.4],popPSTHfsi - std(fineHist(findPVs,:))/sqrt(length(findPVs)),'r','LineWidth',1)
xlim([-0.1 0.3])
title('FSI Population PSTH')
set(gca,'tickdir','out')

subplot(2,1,2)
hold on
plot([-0.2:binSize:0.4],popPSTHmsn,'k','LineWidth',2)
plot([-0.2:binSize:0.4],popPSTHmsn + std(fineHist(findMSNs,:))/sqrt(length(findMSNs)),'k','LineWidth',1)
plot([-0.2:binSize:0.4],popPSTHmsn - std(fineHist(findMSNs,:))/sqrt(length(findMSNs)),'k','LineWidth',1)
xlim([-0.1 0.3])
title('MSN Population PSTH')
set(gca,'tickdir','out')


spikeGraphName = 'PopPSTHMSNPV';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Now lets re-calculate BFs based on binVals
uniqueFreqs = [4000;4594.79342000000;5278.03164300000;6062.86626600000;6964.40450600000;8000;9189.58684000000;10556.0632900000;12125.7325300000;13928.8090100000;16000;18379.1736800000;21112.1265700000;24251.4650600000;27857.6180300000;32000];
newBFs = zeros(size(binValBigStore,2),length(binValBigStore));
bestBF = zeros(length(binValBigStore),1);
for i = 1:length(binValBigStore)
    for j = 1:size(binValBigStore,2)
        [maxVal maxInd] = max(binValBigStore(:,j,i));
        %determine if max is positive
        if maxVal > 0
            newBFs(j,i) = maxInd;
        end
        maxValStore(j) = maxVal;
    end
    %determine best amplitude. 
    [maxValOver maxIndOver] = max(maxValStore);
    if maxVal > 0
        bestBF(i) = newBFs(maxIndOver,i);
    end
end
newBFBack = newBFs;
zeroFind = find(newBFs == 0);
newBFs(zeroFind) = 1;
newBFs = uniqueFreqs(newBFs);
newBFs(zeroFind) = NaN;

bestBFBack = bestBF;
zeroFind = find(bestBF == 0);
bestBF(zeroFind) = 1;
bestBF = uniqueFreqs(bestBF);
bestBF(zeroFind) = NaN;

%% Generate tuning curves that are step functions based on significance.
%how lets look at sigValBigStore
sigValConv = sigValBigStore;
sigValConv(sigValConv <= 0.001) = 4;
sigValConv(sigValConv <= 0.01) = 3;
sigValConv(sigValConv <= 0.05) = 2;
sigValConv(sigValConv <= 1) = 1;
clims = [1 3];


%% Try fitting gaussian?
%first, lets just plot things! Lets just look at what the 70dB band looks
%like
% 
% %plot FSIs
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% plotLim = ceil(sqrt(length(findPVs)));
% figure
% for i = 1:length(findPVs)
%     subplot(plotLim,plotLim,i)
%     plot(squeeze(binValBigStore(:,:,findPVs(i))))
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end
% 
% %plot MSNs
% plotLim = ceil(sqrt(length(findMSNs)));
% figure
% for i = 1:length(findMSNs)
%     subplot(plotLim,plotLim,i)
%     plot(squeeze(binValBigStore(:,:,findMSNs(i))))
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end

%based on this, it looks like one simple amplitude cut might not be the
%best way to do things...lets find the peak amplitude band and average that
%with the one above it. 

%store best amplitude band and extra amplitude band with it. 
ampBand = zeros(length(bigMaster),1);
tarAmpBands = zeros(length(bigMaster),2);
tuneCurveStore = zeros(length(bigMaster),size(binValBigStore,1));
for i = 1:length(bigMaster)
    %find maximum values!
    [C maxInd] =max(binValBigStore(:,:,i));
    %check if all equal
    if all(C == C(1))%case where all equal!
        ampBand(i) = NaN;
        tarAmpBands(i,:) = [NaN,NaN];
    elseif all(C <= 0) %case where all negative
        ampBand(i) = NaN;
        tarAmpBands(i,:) = [NaN,NaN];
    else %if neither of the above cases
        [C2 ampBand(i)] = max(C);
        if ampBand(i) <=4
            tarAmpBands(i,:) = [ampBand(i),ampBand(i) + 1];
        else
            tarAmpBands(i,:) = [ampBand(i) - 1,ampBand(i)];
        end
        %now store data!
        tuneCurveStore(i,:) = mean(binValBigStore(:,tarAmpBands(i,:),i)');
    end 
    
end

% %plot to see how this looks
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% plotLim = ceil(sqrt(length(findPVs)));
% figure
% for i = 1:length(findPVs)
%     subplot(plotLim,plotLim,i)
%     plot(squeeze(tuneCurveStore(findPVs(i),:)))
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end
% 
% %plot MSNs
% plotLim = ceil(sqrt(length(findMSNs)));
% figure
% for i = 1:length(findMSNs)
%     subplot(plotLim,plotLim,i)
%     plot(squeeze(tuneCurveStore(findMSNs(i),:)))
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end

%looks good but isnt the most even. Lets try instead to use sigValConv and
%the threshold amplitude. 

ampBandThresh = zeros(length(bigMaster),1);
for i = 1:length(bigMaster)
    tester = sigValConv(2:end,:,i);
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    if testFind
        ampBandThresh(i) = testFind;
    end
end

%here, i've noticed that a subset of units have significant thresholds but
%a weaker peak based purely on max detection. These units are no good and
%should be eliminated. 

ampBandElim = ampBand - ampBandThresh;
ampBandThresh(ampBandElim<0) = 0;

% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% plotLim = ceil(sqrt(length(findPVs)));
% figure
% for i = 1:length(findPVs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(tuneCurveStore(findPVs(i),:)))
%     if ampBandThresh(findPVs(i)) > 0
%         plot(squeeze(binValBigStore(:,ampBandThresh(findPVs(i)),findPVs(i))))
%         
%     end
%     set(gca,'xtick',[])
%         set(gca,'ytick',[])
% end
% 
% %plot MSNs
% plotLim = ceil(sqrt(length(findMSNs)));
% figure
% for i = 1:length(findMSNs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(tuneCurveStore(findMSNs(i),:)))
%     if ampBandThresh(findMSNs(i)) > 0
%         plot(squeeze(binValBigStore(:,ampBandThresh(findMSNs(i)),findMSNs(i))))
%         
%     end
%     set(gca,'xtick',[])
%         set(gca,'ytick',[])
% end

%Okay, lets fit gaussians to tuneCurveStore

for i = 1:length(tuneCurveStore)
    try
        [fitobject,gof] = fit([1:16]',tuneCurveStore(i,:)','gauss1');
        rSquareStore(i) = gof.rsquare;
    catch
        rSquareStore(i) = 0;
    end
end

% https://www.nature.com/articles/s41467-019-08350-7#Sec9 Chen et al 2019
% use R2 of 0.4. Lets use that here as well. As coincidence detector, use
% big overall significance value (masterData(:,9)) 

findGoodR = find(rSquareStore > 0.4);
findOverSig = find(bigMaster(:,9) == 1);

findGoodRSig = intersect(findGoodR,findOverSig);

% %check this!
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% plotLim = ceil(sqrt(length(findPVs)));
% figure
% for i = 1:length(findPVs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(tuneCurveStore(findPVs(i),:)))
%     if ismember(findPVs(i),findGoodRSig)
%         plot(squeeze(tuneCurveStore(findPVs(i),:)),'r')
%     end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end
% 
% %plot MSNs
% plotLim = ceil(sqrt(length(findMSNs)));
% figure
% for i = 1:length(findMSNs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(tuneCurveStore(findMSNs(i),:)))
%     if ismember(findMSNs(i),findGoodRSig)
%         plot(squeeze(tuneCurveStore(findMSNs(i),:)),'r')
%     end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end

%now lets do something with this??

altPVs = intersect(findPVs,findGoodRSig);
altMSNs = intersect(findMSNs,findGoodRSig);

%try plotting histograms to see how bad things are.


% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% plotLim = ceil(sqrt(length(findPVs)));
% figure
% for i = 1:length(findPVs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(crudeHist(findPVs(i),:)))
%     if ismember(findPVs(i),findGoodRSig)
%         plot(squeeze(crudeHist(findPVs(i),:)),'r')
%     end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end
% 
% %plot MSNs
% plotLim = ceil(sqrt(length(findMSNs)));
% figure
% for i = 1:length(findMSNs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(crudeHist(findMSNs(i),:)))
%     if ismember(findMSNs(i),findGoodRSig)
%         plot(squeeze(crudeHist(findMSNs(i),:)),'r')
%     end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end




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

%% Find some indices
[indPkTr] = functionCellStringFind(masterHeader,'PeakTrough');
[indISI] = functionCellStringFind(masterHeader,'isiCov');
[indBaseFire] = functionCellStringFind(masterHeader,'BaselineRate');
[indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
[indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');

%make vector for width. if neuron responds to everything, that makes 5 *
%16, or 80.
widthHistVect = [0:1:40];


%% Now lets select for units with five significant positive responses. 
signStore = sign(binValBigStore);
signSigStore = sigValBigStore.*signStore;
signSigStore(signSigStore <= 0) = NaN;
sigCut = 0.01;
for i = 1:size(signSigStore,3)
    findSigs(i) = length(find(signSigStore(:,:,i) <= sigCut));
end
% 
% figure
% hist(findSigs,[0:1:20])
% xlim([-0.5 20.5])

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
% 
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% plotLim = ceil(sqrt(length(findPVs)));
% figure
% for i = 1:length(findPVs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(crudeHist(findPVs(i),:)))
%     if ismember(findPVs(i),tarCells)
%         plot(squeeze(crudeHist(findPVs(i),:)),'r')
%     end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end
% 
% %plot MSNs
% plotLim = ceil(sqrt(length(findMSNs)));
% figure
% for i = 1:length(findMSNs)
%     subplot(plotLim,plotLim,i)
%     hold on
%     plot(squeeze(crudeHist(findMSNs(i),:)))
%     if ismember(findMSNs(i),tarCells)
%         plot(squeeze(crudeHist(findMSNs(i),:)),'r')
%     end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
% end

%% NOW LOOK AT LATENCY

latHistVect = [0:0.001:0.1];

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
allZFSI = zHists((tarPVs),:)./max(zHists((tarPVs),:)')';
allZMSN = zHists((tarMSNs),:)./max(zHists((tarMSNs),:)')';

%smooth them
for i = 1:length(tarPVs)
    allZFSI(i,:) = smooth(allZFSI(i,:),3);
end

for i = 1:length(tarMSNs)
    allZMSN(i,:) = smooth(allZMSN(i,:),3);
end

allZFSISTD = std(allZFSI)/sqrt(length(tarPVs));
allZMSNFTD = std(allZMSN)/sqrt(length(tarMSNs));
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI)),'r','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))+allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))-allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN)),'k','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))-allZMSNFTD,'k','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))+allZMSNFTD,'k','LineWidth',1)
set(gca,'TickDir','out');
% xlim([-0.02 0.05])
xlim([0 0.02])
spikeGraphName = '5ValAverageZPlotsMSNrFSIbZOOMNORMALIZED';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
hold on

plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI)),'r','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))+allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))-allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN)),'k','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))-allZMSNFTD,'k','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))+allZMSNFTD,'k','LineWidth',1)

set(gca,'TickDir','out');
% xlim([-0.02 0.05])
xlim([-0.1 0.2])
spikeGraphName = '5ValAverageZPlotsMSNrFSIbNORMALIZED';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

% PLOT WITH STD instead of SEM 

hFig = figure;
hold on
allZFSI = zHists((tarPVs),:)./max(zHists((tarPVs),:)')';
allZMSN = zHists((tarMSNs),:)./max(zHists((tarMSNs),:)')';

%smooth them
for i = 1:length(tarPVs)
    allZFSI(i,:) = smooth(allZFSI(i,:),3);
end

for i = 1:length(tarMSNs)
    allZMSN(i,:) = smooth(allZMSN(i,:),3);
end

allZFSISTD = std(allZFSI);
allZMSNFTD = std(allZMSN);
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI)),'r','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))+allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))-allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN)),'k','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))-allZMSNFTD,'k','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))+allZMSNFTD,'k','LineWidth',1)
set(gca,'TickDir','out');
% xlim([-0.02 0.05])
xlim([0 0.02])
spikeGraphName = '5ValAverageZPlotsMSNrFSIbZOOMNORMALIZEDDSTD';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
hold on

plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI)),'r','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))+allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZFSI)/max(mean(allZFSI))-allZFSISTD,'r','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN)),'k','LineWidth',2)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))-allZMSNFTD,'k','LineWidth',1)
plot([-.2:0.0005:0.4],mean(allZMSN)/max(mean(allZMSN))+allZMSNFTD,'k','LineWidth',1)

set(gca,'TickDir','out');
% xlim([-0.02 0.05])
xlim([-0.1 0.2])
spikeGraphName = '5ValAverageZPlotsMSNrFSIbNORMALIZEDSTD';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%we can also now try and do this for selected frequencies
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

%lets try and do this without spreading in y axis
hFig = figure;
specFind = find(toneSpikes(tarPVs)>0);
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

ranksum(latStore(tarMSNs),latStore(tarPVs))

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

spikeGraphName = '5ValPVMSNSelFreqLatHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%plot as cumulative distribution?
msnLatDist = hist(latStore(tarMSNs),latHistVect);
fsiLatDist = hist(latStore(tarPVs),latHistVect);

msnLatDist = cumsum(msnLatDist);
fsiLatDist = cumsum(fsiLatDist);

latDistNumMSN = msnLatDist(end);
latDistNumFSI = fsiLatDist(end);

msnLatDist = msnLatDist/msnLatDist(end);
fsiLatDist = fsiLatDist/fsiLatDist(end);

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
hold on
plot(latHistVect,msnLatDist,'k')
plot(latHistVect,fsiLatDist,'r')

set(gca,'TickDir','out');
title(['MSN vs FSI Lat Cum Dist-',num2str(latDistNumMSN),'-fSI-',num2str(latDistNumFSI)])

spikeGraphName = '5ValPVMSNSelFreqLatCumDist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%% Now plot heatmaps of target cells
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


%% NOW CALCULATE WIDTHS
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


%% now lets try and plot out width, threshold, with the width of 3 setup.
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
hist(bestBF(tarCells(tarPVs)),[4000;4594.79342000000;5278.03164300000;6062.86626600000;6964.40450600000;8000;9189.58684000000;10556.0632900000;12125.7325300000;13928.8090100000;16000;18379.1736800000;21112.1265700000;24251.4650600000;27857.6180300000;32000],'r')
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
hist(bestBF(tarCells(tarMSNs)),[4000;4594.79342000000;5278.03164300000;6062.86626600000;6964.40450600000;8000;9189.58684000000;10556.0632900000;12125.7325300000;13928.8090100000;16000;18379.1736800000;21112.1265700000;24251.4650600000;27857.6180300000;32000],'k')
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

xlims = [-0.5 17];

trueMax = 10;
for i = 1:5
    
    if finderPV
        testHold = hist(sigWidthPV(i,finderPV),[1:1:16]);
        maxVal = max(testHold);
        if maxVal > trueMax
            trueMax = maxVal;
        end
    end
    
    if finderMSN
        testHold = hist(sigWidthMSN(i,finderMSN),[1:1:16]);
        maxVal = max(testHold);
        if maxVal > trueMax
            trueMax = maxVal;
        end
    end
end

ylims = [0 trueMax + 1;];



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

xlim([0.5 5.5])

spikeGraphName = '5ValPlotMeanWidthChangeAcrossDBs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% now lets plot BF vs width.
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [80 80 1900 1000])
subplot(2,1,1)
plot(bestBF(tarCells(tarPVs)),sigWidthPV(5,:),'k.')
xlabel('BF')
ylabel('Tuning Width')
title('FSI Width vs BF')
set(gca,'TickDir','out');
subplot(2,1,2)
plot(bestBF(tarCells(tarMSNs)),sigWidthMSN(5,:),'k.')
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
 


%% Lets make a nice summary figure

%summary figure will need nine columns, which will be a bit complicated.
%Maybe lets do 10 instead. We want to have 4 rows as well. 

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.04 0.04], [0.04 0.04]);
set(hFig, 'Position', [80 80 1600 1200])

%Row 1, column 2: FSI-MSN separation
subplot(4,5,2)
hold on
plot(peakTrough(:,2)*1000,halfWidth(:,2)*1000,'.','Color',[0.7 0.7 0.7])
plot(peakTrough(findMSNs,2)*1000,halfWidth(findMSNs,2)*1000,'k.')
plot(peakTrough(findPVs,2)*1000,halfWidth(findPVs,2)*1000,'r.')
% plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),halfWidth(bigMaster(:,indCellType) == 2),'g.')
xlabel('Peak Trough (ms)')
ylabel('Half-Width (ms)')
set(gca,'TickDir','out');

%Row 1: column 3-5: some representation of responses of FSI vs MSN. 
load('180622_ML180515D_R17_3300_secondfullTuningFullTuningAnalysis.mat')
tarMSN = 'nt2cluster1';
tarFSI = 'nt16cluster1';

%load frequency specific histograms to create wall of responses. 
histBinVector = s.(tarMSN).HistBinVector;
numDBs = s.SoundData.NumDBs;
numFreqs = s.SoundData.NumFreqs - s.SoundData.WhiteNoise;

%I want 0 before and 0.2 after
targetVect(1) = find(histBinVector< 0,1,'last');
targetVect(2) = find(histBinVector<0.2,1,'last');
% binSize = 0.005;
% reps = 36;
smoothWindWall = 3;

testData = s.(tarMSN);
spacerLength = 10;
histStoreMSN = [];

lengthTrace = length([targetVect(1):targetVect(2)]);
for i = 1:numDBs
    counter = 1;
    for j = 1+s.SoundData.WhiteNoise:s.SoundData.NumFreqs
        histStoreMSN(counter:counter + lengthTrace - 1,i) = smooth(testData.FreqDBHistograms(j,i,[targetVect(1):targetVect(2)]),smoothWindWall);
        counter = counter + lengthTrace;
        histStoreMSN(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
end

maxVal = max(max(histStoreMSN));
subplot(4,5,3)
hold on
for i = 1:numDBs
    plot(histStoreMSN(:,i)/maxVal + i-1,'k')
end
title(maxVal)
xlim([0 length(histStoreMSN)])
set(gca,'TickDir','out');
set(gca,'XTick',[]);
set(gca,'YTickLabel',[]);
ylim([-0.2 numDBs+0.2])

testData = s.(tarFSI);
spacerLength = 10;
histStoreFSI = [];

lengthTrace = length([targetVect(1):targetVect(2)]);
for i = 1:numDBs
    counter = 1;
    for j = 1+s.SoundData.WhiteNoise:s.SoundData.NumFreqs
        histStoreFSI(counter:counter + lengthTrace - 1,i) = smooth(testData.FreqDBHistograms(j,i,[targetVect(1):targetVect(2)]),smoothWindWall);
        counter = counter + lengthTrace;
        histStoreFSI(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
end

maxVal = max(max(histStoreFSI));
subplot(4,5,4)
hold on
for i = 1:numDBs
    plot(histStoreFSI(:,i)/maxVal + i-1,'r')
end
title(maxVal)
xlim([0 length(histStoreFSI)])
set(gca,'TickDir','out');
set(gca,'XTick',[]);
set(gca,'YTickLabel',[]);
ylim([-0.2 numDBs+0.2])



%plot overall histogram MSN
subplot(6,10,9)
meanHist = mean(s.(tarMSN).FrequencyHistograms(1+s.SoundData.WhiteNoise:end,:));
plot(histBinVector,meanHist,'k')
xlim([-0.2 0.4])
set(gca,'TickDir','out');


%plot overall histogram MSN
subplot(6,10,10)
meanHist = mean(s.(tarFSI).FrequencyHistograms(1+s.SoundData.WhiteNoise:end,:));
plot(histBinVector,meanHist,'r')
xlim([-0.2 0.4])
set(gca,'TickDir','out');

%plot waveform MSN
tarWave = s.(tarMSN).AverageWaveForms(:,1);
subplot(12,10,29)
plot(tarWave,'k')
ylim([min(tarWave) max(tarWave)])
xlim([0 length(tarWave)])
set(gca,'YDir','reverse')
set(gca,'XTickLabel',[]);
set(gca,'TickDir','out');

%plot waveform FSI
tarWave = s.(tarFSI).AverageWaveForms(:,1);
subplot(12,10,30)
plot(tarWave,'r')
ylim([min(tarWave) max(tarWave)])
xlim([0 length(tarWave)])
set(gca,'YDir','reverse')
set(gca,'XTickLabel',[]);
set(gca,'TickDir','out');



%Row 2/3, Column 1: plot out pie charts

subplot(4,5,6)
holder = bigMaster(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
% labels(detZero) = [];
% legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('MSNs n=',num2str(length(findMSNs))))

subplot(4,5,11)
holder = bigMaster(findPVs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
% labels(detZero) = [];
% legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('PVs n=',num2str(length(findPVs))))


% Row 2/3 Column 2: pull out threshold

threshMSN = hist(threshStoreMSN,[1:1:5]);
threshFSI = hist(threshStorePV,[1:1:5]);

ylimThresh = max([max(threshMSN),max(threshFSI)]);
subplot(4,5,7)
bar(threshMSN,'k')
xlim([0.5 5.5])
title(strcat('msnThresh u:',num2str(mean(threshStoreMSN))))
ylabel('Number of Cells')
set(gca,'XTickLabel',[]);
ylim([0 ylimThresh])
set(gca,'TickDir','out');

subplot(4,5,12)
bar(threshFSI,'r')
xlim([0.5 5.5])
ylim([0 ylimThresh])
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<30','40','50','60','70'});
title(strcat({'fsiThresh u:',num2str(mean(threshStorePV))},{'pval',num2str(ranksum(threshStorePV,threshStoreMSN))}))
set(gca,'TickDir','out');

%Row 2/3 Column 3: width
widthMSN = hist(widthValsMSN,[1:1:16]);
widthFSI = hist(widthValsPV,[1:1:16]);
ylimWidth = max([max(widthMSN),max(widthFSI)]);

subplot(4,5,8)
bar(widthMSN,'k')
xlim([0.5 16.5])
ylim([0 ylimWidth])
ylabel('Number of Cells')
set(gca,'XTickLabel',[]);
title(strcat('msnwidth u:',num2str(mean(widthValsMSN))))
set(gca,'TickDir','out');


subplot(4,5,13)
bar(widthFSI,'r')
xlim([0.5 16.5])
ylim([0 ylimWidth])
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat({'fsiwidth u:',num2str(mean(widthValsPV))},{'pval',num2str(ranksum(widthValsPV,widthValsMSN))}))
% title('Tuning Width (FSIs)')
set(gca,'TickDir','out');

%Row 2/3 Column 4: latency!
latsMSN = hist(latStore(tarMSNs),latHistVect);
latsFSI = hist(latStore(tarPVs),latHistVect);
latEnd = 0.05;

ylimLats = max([max(latsMSN),max(latsFSI)]);

subplot(4,5,9)
bar(latHistVect,latsMSN,'k')
xlim([latHistVect(1) latEnd])
ylim([0 ylimLats])
set(gca,'TickDir','out');
set(gca,'XTickLabel',[]);
ylabel('Number of Cells')
title(['msnLat u:',num2str(nanmean(latStore(tarMSNs)))])

subplot(4,5,14)
bar(latHistVect,latsFSI,'r')
xlim([latHistVect(1) latEnd])
ylim([0 ylimLats])
set(gca,'TickDir','out');
ylabel('Number of Cells')
title(['fsiLat u:',num2str(nanmean(latStore(tarPVs)))])

%Row 2/3 Column 5, BF
bfsMSN = hist(newBFs(5,tarCells(tarMSNs)),[4000;4594.79342000000;5278.03164300000;6062.86626600000;6964.40450600000;8000;9189.58684000000;10556.0632900000;12125.7325300000;13928.8090100000;16000;18379.1736800000;21112.1265700000;24251.4650600000;27857.6180300000;32000]);
bfsFSI = hist(newBFs(5,tarCells(tarPVs)),[4000;4594.79342000000;5278.03164300000;6062.86626600000;6964.40450600000;8000;9189.58684000000;10556.0632900000;12125.7325300000;13928.8090100000;16000;18379.1736800000;21112.1265700000;24251.4650600000;27857.6180300000;32000]);
bfsYlim = max([max(bfsMSN),max(bfsFSI)]);

subplot(4,5,10)
bar(bfsMSN,'k')
xlim([0.5 16.5])
ylim([0 bfsYlim])
ylabel('Number of Cells')
set(gca,'TickDir','out');
set(gca,'XTickLabel',[]);
title('msnBFs')

subplot(4,5,15)
bar(bfsFSI,'r')
xlim([0.5 16.5])
ylim([0 bfsYlim])
ylabel('Number of Cells')
set(gca,'TickDir','out');
title('fsiBFs')


spikeGraphName = 'NewBaselineOverallFigure';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
print(hFig,spikeGraphName,'-deps','-r0')




%% Supplemental figure: PSTH examples!

%generate for just FSIs.
holder = bigMaster(findPVs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);

negCells = findPVs(det == -2);
mixCells = findPVs(det == -1);
posCells = findPVs(det == 1);


figure
for i = 1:15
subplot(3,5,i)
plot(crudeVect,crudeHist(negCells(i),:))
end
%6 look good. 

figure
for i = 1:15
subplot(3,5,i)
plot(crudeVect,crudeHist(mixCells(i),:))
end
%14  work

figure
for i = 1:15
subplot(3,5,i)
plot(crudeVect,crudeHist(posCells(i),:))
end
%1, 24, 25 look good. 

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.04 0.04], [0.04 0.04]);
set(hFig, 'Position', [80 80 1600 1200])

subplot(1,3,1)
plot(crudeVect,smooth(crudeHist(findPVs(1),:),5))
xlim([-0.15 0.35])
subplot(1,3,2)
plot(crudeVect,smooth(crudeHist(findPVs(81),:),5))
xlim([-0.15 0.35])
subplot(1,3,3)
plot(crudeVect,smooth(crudeHist(tarCells(tarPVs(10)),:),5))
xlim([-0.15 0.35])


spikeGraphName = 'ExamplePSTHs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
print(hFig,spikeGraphName,'-deps','-r0')


%generate from just MSNs. 

holder = bigMaster(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);

nonCells = findMSNs(det == 0);
negCells = findMSNs(det == -2);
mixCells = findMSNs(det == -1);
posCells = findMSNs(det == 1);


figure
for i = 1:25
subplot(5,5,i)
plot(crudeVect,crudeHist(negCells(i),:))
end
%21, 23, 24 look good. 

figure
for i = 1:25
subplot(5,5,i)
plot(crudeVect,crudeHist(mixCells(i),:))
end
%11 and 12 both work

figure
for i = 1:25
subplot(5,5,i)
plot(crudeVect,crudeHist(posCells(i),:))
end
%11, 19, 25 look good. 



hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.04 0.04], [0.04 0.04]);
set(hFig, 'Position', [80 80 1600 1200])
subplot(2,2,1)
plot(crudeVect,smooth(crudeHist(nonCells(9),:),5))
xlim([-0.15 0.35])
ylim([0 2.5])
set(gca,'TickDir','out')
subplot(2,2,2) %positive
plot(crudeVect,smooth(crudeHist(posCells(11),:),5))

xlim([-0.15 0.35])
ylim([0 2.5])
set(gca,'TickDir','out')
subplot(2,2,3) %negative
plot(crudeVect,smooth(crudeHist(negCells(21),:),5))

xlim([-0.15 0.35])
ylim([0 2.5])
set(gca,'TickDir','out')
subplot(2,2,4) %mixed
plot(crudeVect,smooth(crudeHist(mixCells(11),:),5))
xlim([-0.15 0.35])
ylim([0 2.5])
set(gca,'TickDir','out')


spikeGraphName = 'ExamplePSTHsMSN';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
print(hFig,spikeGraphName,'-deps','-r0')



%% Generate population tuning curves
widthPer = 0.5;

for bigInd = 1:5
    tmpDB = bigInd;

    bigStoreMSN = NaN(length(tarMSNs),41);
    bigStoreFSI = NaN(length(tarPVs),41);

    for i = 1:length(tarMSNs)
        %pull curves
        tmp1 = binValBigStore(:,tmpDB,tarCells(tarMSNs(i)));
        %determine relative max.
        [maxVal1 maxInd1] = max(tmp1);
        %save max index. 
        bfSaveMSN(i,bigInd) = maxInd1;
        %divide by maximum value
        if maxVal1 > 0
            bigStoreMSN(i,21-maxInd1+1:21-maxInd1+length(tmp1)) = tmp1/maxVal1;
        end
    end
    
    for i = 1:length(tarPVs)
        %pull curves
        tmp1 = binValBigStore(:,tmpDB,tarCells(tarPVs(i)));
        %determine relative max.
        [maxVal1 maxInd1] = max(tmp1);
        %save max index. 
        bfSaveFSI(i,bigInd) = maxInd1;
        %divide by maximum value
        if maxVal1 > 0
            bigStoreFSI(i,21-maxInd1+1:21-maxInd1+length(tmp1)) = tmp1/maxVal1;
        end
    end

    hFig = figure;
    subplot(2,1,1)
    plot(nanmean(bigStoreMSN),'k')
    hold on
    plot(nanmean(bigStoreFSI),'r')
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:21) <= widthPer,1,'last');
    widthMSN = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreFSI);
    find1 = find(tester(1:21) <= widthPer,1,'last');
    widthFSI = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tmpDB),'widthMSN:',num2str(widthMSN),'widthFSI:',num2str(widthFSI)])
    set(gca,'TickDir','out')
    xlim([1 41])
    
    
    subplot(2,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreFSI);
    test = sum(tester);
    plot(test,'g')
    xlim([1 41])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    
    spikeGraphName = strcat(['AverageTuningAlignedToBaseAtDB-',num2str(tmpDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end


%% Now apply this restricted analysis to FSIs. May not change much. 

for bigInd = 1:5
    tmpDB = bigInd;

    bigStoreMSN = NaN(length(tarMSNs),41);
    bigStoreFSI = NaN(length(tarPVs),41);

    for i = 1:length(tarMSNs)
        if sigStore5Val(bigInd,(tarMSNs(i))) == 1
            %pull curves
            tmp1 = binValBigStore(:,tmpDB,tarCells(tarMSNs(i)));
            %determine relative max.
            [maxVal1 maxInd1] = max(tmp1);
            %save max index. 
            bfSaveMSN(i,bigInd) = maxInd1;
            %divide by maximum value
            if maxVal1 > 0
                bigStoreMSN(i,21-maxInd1+1:21-maxInd1+length(tmp1)) = tmp1/maxVal1;
            end
        end
    end
    
    for i = 1:length(tarPVs)
        if sigStore5Val(bigInd,(tarPVs(i))) == 1
            %pull curves
            tmp1 = binValBigStore(:,tmpDB,tarCells(tarPVs(i)));
            %determine relative max.
            [maxVal1 maxInd1] = max(tmp1);
            %save max index. 
            bfSaveFSI(i,bigInd) = maxInd1;
            %divide by maximum value
            if maxVal1 > 0
                bigStoreFSI(i,21-maxInd1+1:21-maxInd1+length(tmp1)) = tmp1/maxVal1;
            end
        end
    end

    hFig = figure;
    subplot(2,1,1)
    hold on
    plot(nanmean(bigStoreMSN),'k','LineWidth',2)
    plot(nanmean(bigStoreMSN) + nanstd(bigStoreMSN)/sqrt(length(tarMSNs)),'k','LineWidth',1)
    plot(nanmean(bigStoreMSN) - nanstd(bigStoreMSN)/sqrt(length(tarMSNs)),'k','LineWidth',1)
    plot(nanmean(bigStoreFSI),'r')
    plot(nanmean(bigStoreFSI) + nanstd(bigStoreFSI)/sqrt(length(tarPVs)),'r','LineWidth',1)
    plot(nanmean(bigStoreFSI) - nanstd(bigStoreFSI)/sqrt(length(tarPVs)),'r','LineWidth',1)
    tester = nanmean(bigStoreMSN);
    find1 = find(tester(1:21) <= widthPer,1,'last');
    widthMSN = find(tester(find1+1:end) <= widthPer,1,'first');
    tester = nanmean(bigStoreFSI);
    find1 = find(tester(1:21) <= widthPer,1,'last');
    widthFSI = find(tester(find1+1:end) <= widthPer,1,'first');
    %find width at specified percentage
    title(['AvAlignedToBaseAtDB-',num2str(tmpDB),'widthMSN:',num2str(widthMSN),'widthFSI:',num2str(widthFSI)])
    set(gca,'TickDir','out')
    xlim([1 41])
    
    
    subplot(2,1,2)
    hold on
    tester = ~isnan(bigStoreMSN);
    test = sum(tester);
    plot(test,'k')
    tester = ~isnan(bigStoreFSI);
    test = sum(tester);
    plot(test,'g')
    xlim([1 41])
    title('Number of Data Points')
    set(gca,'TickDir','out')
    
    spikeGraphName = strcat(['AverageRestrictedTuningAlignedToBaseAtDB-',num2str(tmpDB)]);
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

%% Generate population tuning curve?



popNormTuningMSN
popNormTuningFSI










