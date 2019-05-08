
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
bigpspk = [];
bigpx = [];
bigpxspk = [];

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
        
        %now extract non-linearities
        bigpspk(bigMasterInd + j - 1) = s.NonLinear.pspk{j};
        bigpx(bigMasterInd + j - 1,:) = s.NonLinear.px{j};
        bigpxspk(bigMasterInd + j - 1,:) = s.NonLinear.pxspk{j};

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
        crudeHist(bigMasterInd + j - 1,:) = s.(desigName{j}).AllHistograms;
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


%% FIND SIGNIFICANT STAs

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
dmrMSN = intersect(dmrMSN,selMSNs);
dmrPV = findPVs;
dmrPV(ismember(dmrPV,fewSpikes)) = [];
dmrPV = intersect(dmrPV,sigUnitPrct);
dmrPV = intersect(dmrPV,selPVs);


%% CALCULATE CHANGE FOR FIRING RATE AND PLI
%plot out spiking rate for DMR and PLI

%calculate phase locking index
% PLI = (max(bigSTAstore') - min(bigSTAstore'))./(dmrSpikeNum^2/(10*60)*sqrt(8));
PLI = (max(bigSTAstore') - min(bigSTAstore'))./(dmrSpikeNum*38.8520);


hFig = figure;
set(hFig, 'Position', [10 80 1200 800])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
histVectSpikeNum = [0:1:40];
subplot(2,3,1)
hist(dmrSpikeNum(dmrPV)/(10*60),histVectSpikeNum)
title(['Hist Mean FR FSIs n=',num2str(length(dmrPV))])
subplot(2,3,4)
hist(dmrSpikeNum(dmrMSN)/(10*60),histVectSpikeNum)
title(['Hist Mean FR MSNs n=',num2str(length(dmrMSN))])
%plot out spike rate minus baseline rate from tones
histVectSpikeNum = [-1:0.1:1];
subplot(2,3,2)
hist((dmrSpikeNum(dmrPV)/(10*60)- bigMaster(dmrPV,8)')./(dmrSpikeNum(dmrPV)/(10*60)+ bigMaster(dmrPV,8)'),histVectSpikeNum)
title('Hist Mod Index FSIs')
subplot(2,3,5)
hist((dmrSpikeNum(dmrMSN)/(10*60)- bigMaster(dmrMSN,8)')./(dmrSpikeNum(dmrMSN)/(10*60)+ bigMaster(dmrMSN,8)'),histVectSpikeNum)
title('Hist Mod Index MSNs')

%now plot out PLI
histVectSpikeNum = [0:0.05:1];
subplot(2,3,3)
hist(PLI(dmrPV),histVectSpikeNum)
title('PLI FSI')
subplot(2,3,6)
hist(PLI(dmrMSN),histVectSpikeNum)
title('PLI MSN')

spikeGraphName = 'STA FR Change FR PLI';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% plot out STAs

map = [0 0 1];
for i = 2:21
    if i < 11
        map(i,:) = map(i-1,:);
        map(i,1:2) = map(i,1:2) + 0.1;
    elseif i == 11
        map(i,:) = [1 1 1];
    elseif i > 11
        map(i,:) = map(i-1,:);
        map(i,2:3) = map(i,2:3) - 0.1;
    end
end

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.01 0.03], [0.05 0.05]);
axisVal = ceil(sqrt(length(dmrPV)));
for i = 1:length(dmrPV)
    subplot(axisVal,axisVal,i)
    limVal = max(abs(bigSTASigstore(dmrPV(i),:)));
    clims = [-1.05 * limVal,1.05 * limVal];
    imagesc(reshape(bigSTASigstore(dmrPV(i),:),length(faxis),[]),clims)
    colormap(map)
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    title(num2str(dmrSpikeNum(dmrPV(i))))
end

spikeGraphName = 'fsiSTAs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.01 0.03], [0.05 0.05]);
axisVal = ceil(sqrt(length(dmrMSN)));
for i = 1:length(dmrMSN)
    subplot(axisVal,axisVal,i)
    limVal = max(abs(bigSTASigstore(dmrMSN(i),:)));
    clims = [-1.05 * limVal,1.05 * limVal];
    imagesc(reshape(bigSTASigstore(dmrMSN(i),:),length(faxis),[]),clims)
    colormap(map)
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    title(num2str(dmrSpikeNum(dmrMSN(i))))
end

spikeGraphName = 'msnSTAs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


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

%% to determine BF, lets compress along time axis. 
smoothWind = 7;
% smoothWind = 1;
for i = 1:length(bigDBStore)
    compTimingPos(:,i) = sum(posSTA(:,:,i));
    compTimingNeg(:,i) =  sum(negSTA(:,:,i));
    compWidthPos(:,i) = smooth(sum(posSTA(:,:,i)'),smoothWind);
    compWidthNeg(:,i) = smooth(sum(negSTA(:,:,i)'),smoothWind);
    
end



%Lets examine BFs. Restrict search to below 32 kHz
limiter = [1 31];
[maxVal staBFPos] = max(compWidthPos(limiter(1):limiter(2),:));
staBFPos = faxis(staBFPos);

[maxVal staBFNeg] = min(compWidthNeg(limiter(1):limiter(2),:));
staBFNeg = faxis(staBFNeg);

%BF store from master data is sketchy. Pull directly from final amplitude
%of binDiff. 

[maxVal toneBF] = max(squeeze(sum(squeeze(binValBigStore(:,end-2:end,:)),2)));
uniqueFreqs = s.SoundData.UniqueFrequencies;
toneBF = uniqueFreqs(toneBF);

hFig = figure;
set(hFig, 'Position', [10 80 1500 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.07 0.07], [0.05 0.05]);
subplot(2,3,1)
hold on
loglog(staBFPos(dmrPV),toneBF(dmrPV),'r.')
plot([4000 64000],[4000 64000],'k')
set(gca,'TickDir','out')
xlabel('STA BF (Hz)')
ylabel('Pure Tone BF (Hz)')
title('FSI Positive BF Comparison')

subplot(2,3,4)
hold on
loglog(staBFPos(dmrMSN),toneBF(dmrMSN),'r.')
plot([4000 64000],[4000 64000],'k')
set(gca,'TickDir','out')
xlabel('STA BF (Hz)')
ylabel('Pure Tone BF (Hz)')
title('MSN Positive BF Comparison')

subplot(2,3,2)
hold on
loglog(staBFNeg(dmrPV),toneBF(dmrPV),'r.')
plot([4000 64000],[4000 64000],'k')
set(gca,'TickDir','out')
xlabel('STA BF (Hz)')
ylabel('Pure Tone BF (Hz)')
title('FSI Negative BF Comparison')

subplot(2,3,5)
hold on
loglog(staBFNeg(dmrMSN),toneBF(dmrMSN),'r.')
plot([4000 64000],[4000 64000],'k')
set(gca,'TickDir','out')
xlabel('STA BF (Hz)')
ylabel('Pure Tone BF (Hz)')
title('MSN Negative BF Comparison')


subplot(2,3,3)
hold on
loglog(staBFNeg(dmrPV),staBFPos(dmrPV),'r.')
plot([4000 64000],[4000 64000],'k')
set(gca,'TickDir','out')
xlabel('STA BF (Neg) (Hz)')
ylabel('STA BF (Pos) (Hz)')
title('FSI STA BF Comparison')

subplot(2,3,6)
hold on
loglog(staBFNeg(dmrMSN),staBFPos(dmrMSN),'r.')
plot([4000 64000],[4000 64000],'k')
set(gca,'TickDir','out')
xlabel('STA BF (Neg) (Hz)')
ylabel('STA BF (Pos) (Hz)')
title('MSN STA BF Comparison')

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

%% Look at widths 
dmrWidthPos = [];
dmrBFPos = [];
dmrWidthNeg = [];
dmrBFNeg = [];
widthPer = 0.5;
for i = 1:length(bigDBStore)
    %pull peak modulation timing. 
    [maxVal dmrTimePos(i)] = max(compTimingPos(:,i));
    [minVal dmrTimeNeg(i)] = min(compTimingPos(:,i));
    %pull widths from STA
    [dmrWidthPos(i,1:4),dmrBFPos(i),maxVal,cutVal] = functionHeightBasedTuningWidth(compWidthPos(:,i),widthPer,1);
    [dmrWidthNeg(i,1:4),dmrBFNeg(i),maxVal,cutVal] = functionHeightBasedTuningWidth(-1*compWidthNeg(:,i),widthPer,1);
    %pull width from tuning
    toneWidth = binValBigStore(:,:,i);
    toneWidth = mean(toneWidth(:,end-2:end)');
    %extract width!
    [toneWidthPos(i,1:4),toneBFPos(i),maxVal,cutVal] = functionHeightBasedTuningWidth(toneWidth,widthPer,1);
end

dmrWidthPos (:,[2,4]) = [];
dmrWidthNeg (:,[2,4]) = [];
toneWidthPos (:,[2,4]) = [];


%first, lets just make plots of tuning curves of pure tones vs DMR.
%start with PVs
axisDef = ceil(sqrt(length(dmrPV)));
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
for i = 1:length(dmrPV)
    %pull tuning curve average from top three amplitudes
    toneWidth = binValBigStore(:,:,dmrPV(i));
    toneWidth = mean(toneWidth(:,end-2:end)');
    %make subplot
    subplot(axisDef,axisDef,i)
    hold on
    %plot tone information. General curve
    plot(toneFreqs,toneWidth/max(toneWidth),'k')
    %plot out peak value
    plot(toneFreqs(toneBFPos(dmrPV(i))),toneWidth(toneBFPos(dmrPV(i)))/max(toneWidth),'r*')
    %plot out width
    if ~isnan(toneWidthPos(dmrPV(i),1))
        width1 = toneBFPos(dmrPV(i)) - toneWidthPos(dmrPV(i),1);
        width1 = interp1([1:16],toneFreqs,width1);
    else
        width1 = toneFreqs(1);
    end
    
    if ~isnan(toneWidthPos(dmrPV(i),2))
        width2 = toneBFPos(dmrPV(i)) + toneWidthPos(dmrPV(i),2);
        width2 = interp1([1:16],toneFreqs,width2);
    else
        width2 = toneFreqs(end);
    end
    plot([width1 width2],[0.5 0.5],'r','LineWidth',2)
    
    %now plot STA info!    
    plot(faxis,compWidthPos(:,dmrPV(i))/max(compWidthPos(:,dmrPV(i))),'r')
    %plot out BF
    plot(faxis(dmrBFPos(dmrPV(i))),1,'g*')
    %plot out width!
    if ~isnan(dmrWidthPos(dmrPV(i),1))
        width1 = dmrBFPos(dmrPV(i)) - dmrWidthPos(dmrPV(i),1);
        width1 = interp1([1:40],faxis,width1);
    else
        width1 = faxis(1);
    end
    
    if ~isnan(dmrWidthPos(dmrPV(i),2))
        width2 = dmrBFPos(dmrPV(i)) + dmrWidthPos(dmrPV(i),2);
        width2 = interp1([1:40],faxis,width2);
    else
        width2 = faxis(end);
    end
    plot([width1 width2],[0.4 0.4],'g','LineWidth',2)
    %plot negative!
    plot(faxis,-1*compWidthNeg(:,dmrPV(i))/min(compWidthNeg(:,dmrPV(i))),'b')
    %plot max!
    plot(faxis(dmrBFNeg(dmrPV(i))),-1,'g*')
    %plot width!
    if ~isnan(dmrWidthNeg(dmrPV(i),1))
        width1 = dmrBFNeg(dmrPV(i)) - dmrWidthNeg(dmrPV(i),1);
        width1 = interp1([1:40],faxis,width1);
    else
        width1 = faxis(1);
    end
    
    if ~isnan(dmrWidthNeg(dmrPV(i),2))
        width2 = dmrBFNeg(dmrPV(i)) + dmrWidthNeg(dmrPV(i),2);
        width2 = interp1([1:40],faxis,width2);
    else
        width2 = faxis(end);
    end
    plot([width1 width2],[-0.4 -0.4],'g','LineWidth',2)
    
    xlim([faxis(1) faxis(end)])
    set(gca, 'XScale', 'log')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

spikeGraphName = 'FSIFreqWidths';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%Then MSNs
axisDef = ceil(sqrt(length(dmrMSN)));
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
for i = 1:length(dmrMSN)
    %pull tuning curve average from top three amplitudes
    toneWidth = binValBigStore(:,:,dmrMSN(i));
    toneWidth = mean(toneWidth(:,end-2:end)');
    %make subplot
    subplot(axisDef,axisDef,i)
    hold on
    %plot tone information. General curve
    plot(toneFreqs,toneWidth/max(toneWidth),'k')
    %plot out peak value
    if ~isnan(toneBFPos(dmrMSN(i)))
        plot(toneFreqs(round(toneBFPos(dmrMSN(i)))),toneWidth(round(toneBFPos(dmrMSN(i))))/max(toneWidth),'r*')
        %plot out width
        if ~isnan(toneWidthPos(dmrMSN(i),1))
            width1 = toneBFPos(dmrMSN(i)) - toneWidthPos(dmrMSN(i),1);
            width1 = interp1([1:16],toneFreqs,width1);
        else
            width1 = toneFreqs(1);
        end

        if ~isnan(toneWidthPos(dmrMSN(i),2))
            width2 = toneBFPos(dmrMSN(i)) + toneWidthPos(dmrMSN(i),2);
            width2 = interp1([1:16],toneFreqs,width2);
        else
            width2 = toneFreqs(end);
        end
        plot([width1 width2],[0.5 0.5],'r','LineWidth',2)
    end
    
    
    %now plot STA info!    
    plot(faxis,compWidthPos(:,dmrMSN(i))/max(compWidthPos(:,dmrMSN(i))),'r')
    %plot out BF
    plot(faxis(dmrBFPos(dmrMSN(i))),1,'g*')
    %plot out width!
    if ~isnan(dmrWidthPos(dmrMSN(i),1))
        width1 = dmrBFPos(dmrMSN(i)) - dmrWidthPos(dmrMSN(i),1);
        width1 = interp1([1:40],faxis,width1);
    else
        width1 = faxis(1);
    end
    
    if ~isnan(dmrWidthPos(dmrMSN(i),2))
        width2 = dmrBFPos(dmrMSN(i)) + dmrWidthPos(dmrMSN(i),2);
        width2 = interp1([1:40],faxis,width2);
    else
        width2 = faxis(end);
    end
    plot([width1 width2],[0.4 0.4],'g','LineWidth',2)
    %plot negative!
    plot(faxis,-1*compWidthNeg(:,dmrMSN(i))/min(compWidthNeg(:,dmrMSN(i))),'b')
    %plot max!
    plot(faxis(dmrBFNeg(dmrMSN(i))),-1,'g*')
    %plot width!
    if ~isnan(dmrWidthNeg(dmrMSN(i),1))
        width1 = dmrBFNeg(dmrMSN(i)) - dmrWidthNeg(dmrMSN(i),1);
        width1 = interp1([1:40],faxis,width1);
    else
        width1 = faxis(1);
    end
    
    if ~isnan(dmrWidthNeg(dmrMSN(i),2))
        width2 = dmrBFNeg(dmrMSN(i)) + dmrWidthNeg(dmrMSN(i),2);
        width2 = interp1([1:40],faxis,width2);
    else
        width2 = faxis(end);
    end
    plot([width1 width2],[-0.4 -0.4],'g','LineWidth',2)
    
    xlim([faxis(1) faxis(end)])
    set(gca, 'XScale', 'log')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

spikeGraphName = 'MSNFreqWidths';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets calculate widths based on half height. 
heightBench = 0.5;
%From pure tone:

%calculate widths with height based system
bigWidthHeightStore = [];
for i = 1:length(bigMaster)
    testValues = squeeze(binValBigStore(:,:,i));
    for j = 1:size(testValues,2)
        [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues(:,j),heightBench,3);
        bigWidthHeightStore(j,:,i) = widthVals;
        bigWidthMaxPosStore(j,i) = maxPos;
        bigWidthMaxValStore(j,i) = maxVal;
        bigWidthCutVal(j,i) = cutVal;
    end
end

%now lets just pull out the outward search, since this looks to be far more
%accurate
bigWidthSelWidth = bigWidthHeightStore(:,[1,3],:);
tarCells = [1:length(bigMaster)];
tarPVs = findPVs;
tarMSNs = findMSNs;
%now lets pull the height based width values, and the centers of these
%"responses"
tarCellHeightWidth = bigWidthSelWidth(:,:,tarCells);
tarCellMaxPos = bigWidthMaxPosStore(:,tarCells);


sigWidthStore = [];
sideWarnLow = [];
sideWarnHi = [];
for i = 1:length(tarCells)
    for j = 1:size(tarCellHeightWidth,1)
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


bigStoreSTAPosWidth = [];
for i = 1:length(bigMaster)
    testValues = squeeze(compWidthPos(:,i));
    [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues,heightBench,3);
    bigStoreSTAPosWidth(:,i) = widthVals;
    bigStoreSTAPosWidthPeakFreq(i) = maxPos;
    bigStoreSTAPosWidthPeakVal(i) = maxVal;
    bigStoreSTAPosWidthCutVal(i) = cutVal;
end

%now lets just pull out the outward search, since this looks to be far more
%accurate
bigStoreSTAPosWidthSel = bigStoreSTAPosWidth([1,3],:);

tarCells = [1:length(bigMaster)];
tarPVs = findPVs;
tarMSNs = findMSNs;
%now lets pull the height based width values, and the centers of these
%"responses"
tarCellHeightWidth = bigStoreSTAPosWidthSel(:,tarCells);
tarCellMaxPos = bigStoreSTAPosWidthPeakFreq(tarCells);


sigWidthStoreSTAPos = [];
sideWarnLowSTAPos = [];
sideWarnHiSTAPos = [];
for i = 1:length(tarCells)
    %now check for width values. If NaN on either side, determine
    %edge.
    %start with lower edge. 
    if isnan(tarCellHeightWidth(1,i))
        lowWidth = tarCellMaxPos(i);
        sideWarnLowSTAPos(i) = 1;
    else
        lowWidth = tarCellHeightWidth(1,i);
    end
    if isnan(tarCellHeightWidth(2,i))
        hiWidth = length(faxis) - tarCellMaxPos(i);
        sideWarnHiSTAPos(i) = 1;
    else
        hiWidth = tarCellHeightWidth(2,i);
    end
    sigWidthStoreSTAPos(i) = hiWidth + lowWidth;
end

%compare max amplitude pure tones to STAs

ampTar = 5;

%octaves per step??
unitsTone = 0.2;
unitsSTA = 4/39;

widthCompFSI = [sigWidthStore(ampTar,dmrPV)*unitsTone;sigWidthStoreSTAPos(dmrPV)*unitsSTA];
widthCompMSN = [sigWidthStore(ampTar,dmrMSN)*unitsTone;sigWidthStoreSTAPos(dmrMSN)*unitsSTA];

hFig = figure;
set(hFig, 'Position', [10 80 1000 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);

subplot(2,1,1)
hold on
for i = 1:length(dmrPV)
    plot([1 2],[widthCompFSI(1,i) widthCompFSI(2,i)],'Color',[0.7 0.7 0.7])
end
errorbar([1 2],[nanmean(widthCompFSI(1,:)) nanmean(widthCompFSI(2,:))],[nanstd(widthCompFSI(1,:))/sqrt(length(dmrPV)) nanstd(widthCompFSI(2,:))/sqrt(length(dmrPV))],'k','LineWidth',2)
xlim([0.5 2.5])
set(gca,'TickDir','out');
title('Change in FSI HalfWidth')
subplot(2,1,2)
hold on
for i = 1:length(dmrMSN)
    plot([1 2],[widthCompMSN(1,i) widthCompMSN(2,i)],'Color',[0.7 0.7 0.7])
end
errorbar([1 2],[nanmean(widthCompMSN(1,:)) nanmean(widthCompMSN(2,:))],[nanstd(widthCompMSN(1,:))/sqrt(length(dmrMSN)) nanstd(widthCompMSN(2,:))/sqrt(length(dmrMSN))],'k','LineWidth',2)
xlim([0.5 2.5])
set(gca,'TickDir','out');
title('Change in MSN HalfWidth')

spikeGraphName = 'HalfWidthSTAvsPureTone';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


% %Now lets go to 0.2 of height
% %now lets calculate widths based on half height. 
% heightBench = 0.2;
% %From pure tone:
% 
% %calculate widths with height based system
% bigWidthHeightStore = [];
% for i = 1:length(bigMaster)
%     testValues = squeeze(binValBigStore(:,:,i));
%     for j = 1:size(testValues,2)
%         [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues(:,j),heightBench,3);
%         bigWidthHeightStore(j,:,i) = widthVals;
%         bigWidthMaxPosStore(j,i) = maxPos;
%         bigWidthMaxValStore(j,i) = maxVal;
%         bigWidthCutVal(j,i) = cutVal;
%     end
% end
% 
% %now lets just pull out the outward search, since this looks to be far more
% %accurate
% bigWidthSelWidth = bigWidthHeightStore(:,[1,3],:);
% tarCells = [1:length(bigMaster)];
% tarPVs = findPVs;
% tarMSNs = findMSNs;
% %now lets pull the height based width values, and the centers of these
% %"responses"
% tarCellHeightWidth = bigWidthSelWidth(:,:,tarCells);
% tarCellMaxPos = bigWidthMaxPosStore(:,tarCells);
% 
% 
% sigWidthStore = [];
% sideWarnLow = [];
% sideWarnHi = [];
% for i = 1:length(tarCells)
%     for j = 1:size(tarCellHeightWidth,1)
%         %now check for width values. If NaN on either side, determine
%         %edge.
%         %start with lower edge. 
%         if isnan(tarCellHeightWidth(j,1,i))
%             lowWidth = tarCellMaxPos(j,i);
%             sideWarnLow(j,i) = 1;
%         else
%             lowWidth = tarCellHeightWidth(j,1,i);
%         end
%         if isnan(tarCellHeightWidth(j,2,i))
%             hiWidth = 16 - tarCellMaxPos(j,i);
%             sideWarnHi(j,i) = 1;
%         else
%             hiWidth = tarCellHeightWidth(j,2,i);
%         end
%         sigWidthStore(j,i) = hiWidth + lowWidth;
%     end
% end
% 
% 
% bigStoreSTAPosWidth = [];
% for i = 1:length(bigMaster)
%     testValues = squeeze(compWidthPos(:,i));
%     [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues,heightBench,3);
%     bigStoreSTAPosWidth(:,i) = widthVals;
%     bigStoreSTAPosWidthPeakFreq(i) = maxPos;
%     bigStoreSTAPosWidthPeakVal(i) = maxVal;
%     bigStoreSTAPosWidthCutVal(i) = cutVal;
% end
% 
% %now lets just pull out the outward search, since this looks to be far more
% %accurate
% bigStoreSTAPosWidthSel = bigStoreSTAPosWidth([1,3],:);
% 
% tarCells = [1:length(bigMaster)];
% tarPVs = findPVs;
% tarMSNs = findMSNs;
% %now lets pull the height based width values, and the centers of these
% %"responses"
% tarCellHeightWidth = bigStoreSTAPosWidthSel(:,tarCells);
% tarCellMaxPos = bigStoreSTAPosWidthPeakFreq(tarCells);
% 
% 
% sigWidthStoreSTAPos = [];
% sideWarnLowSTAPos = [];
% sideWarnHiSTAPos = [];
% for i = 1:length(tarCells)
%     %now check for width values. If NaN on either side, determine
%     %edge.
%     %start with lower edge. 
%     if isnan(tarCellHeightWidth(1,i))
%         lowWidth = tarCellMaxPos(i);
%         sideWarnLowSTAPos(i) = 1;
%     else
%         lowWidth = tarCellHeightWidth(1,i);
%     end
%     if isnan(tarCellHeightWidth(2,i))
%         hiWidth = length(faxis) - tarCellMaxPos(i);
%         sideWarnHiSTAPos(i) = 1;
%     else
%         hiWidth = tarCellHeightWidth(2,i);
%     end
%     sigWidthStoreSTAPos(i) = hiWidth + lowWidth;
% end
% 
% %compare max amplitude pure tones to STAs
% 
% ampTar = 5;
% 
% %octaves per step??
% unitsTone = 0.2;
% unitsSTA = 4/39;
% 
% widthCompFSI = [sigWidthStore(ampTar,dmrPV)*unitsTone;sigWidthStoreSTAPos(dmrPV)*unitsSTA];
% widthCompMSN = [sigWidthStore(ampTar,dmrMSN)*unitsTone;sigWidthStoreSTAPos(dmrMSN)*unitsSTA];
% 
% hFig = figure;
% set(hFig, 'Position', [10 80 1000 1000])
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
% 
% subplot(2,1,1)
% hold on
% for i = 1:length(dmrPV)
%     plot([1 2],[widthCompFSI(1,i) widthCompFSI(2,i)],'Color',[0.7 0.7 0.7])
% end
% errorbar([1 2],[nanmean(widthCompFSI(1,:)) nanmean(widthCompFSI(2,:))],[nanstd(widthCompFSI(1,:))/sqrt(length(dmrPV)) nanstd(widthCompFSI(2,:))/sqrt(length(dmrPV))],'k','LineWidth',2)
% xlim([0.5 2.5])
% set(gca,'TickDir','out');
% title('Change in FSI 20% Width')
% subplot(2,1,2)
% hold on
% for i = 1:length(dmrMSN)
%     plot([1 2],[widthCompMSN(1,i) widthCompMSN(2,i)],'Color',[0.7 0.7 0.7])
% end
% errorbar([1 2],[nanmean(widthCompMSN(1,:)) nanmean(widthCompMSN(2,:))],[nanstd(widthCompMSN(1,:))/sqrt(length(dmrMSN)) nanstd(widthCompMSN(2,:))/sqrt(length(dmrMSN))],'k','LineWidth',2)
% xlim([0.5 2.5])
% set(gca,'TickDir','out');
% title('Change in MSN 20% Width')
% 
% spikeGraphName = 'Width20%STAvsPureTone';
% savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')






%% now lets plot out temporal information. 

latFineHist(:,1) = 0;
latFineHist(:,end) = 0;
for i = 1:length(bigpspk)
    %extract std and mean from baseline
    meanVal = mean(latFineHist(tarCells(i),1:401));
    stdVal = std(latFineHist(tarCells(i),1:401));
    allZHists(i,:) = (latFineHist(tarCells(i),:) - meanVal)/stdVal;
end

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


%lets try and pull peak positive mod and peak negative mod. 
[maxVal peakTimePos] = max(compTimingPos);
peakTimePos = (100-peakTimePos)*timeStep;

[maxVal peakTimeNeg] = min(compTimingNeg);
peakTimeNeg = (100-peakTimeNeg)*timeStep;

hFig = figure;
set(hFig, 'Position', [10 80 1000 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(2,2,1)
hist(peakTimePos(dmrPV),[0:0.005:0.1])
xlabel('Time (sec)')
ylabel('Number of Units')
set(gca,'TickDir','out');
title('FSI Latency (Positive Response)')
subplot(2,2,3)
hist(peakTimePos(dmrMSN),[0:0.005:0.1])
xlabel('Time (sec)')
ylabel('Number of Units')
set(gca,'TickDir','out');
title('MSN Latency (Positive Response)')

subplot(2,2,2)
hist(peakTimeNeg(dmrPV),[0:0.005:0.1])
xlabel('Time (sec)')
ylabel('Number of Units')
set(gca,'TickDir','out');
title('FSI Latency (Negative Response)')
subplot(2,2,4)
hist(peakTimeNeg(dmrMSN),[0:0.005:0.1])
xlabel('Time (sec)')
ylabel('Number of Units')
set(gca,'TickDir','out');
title('MSN Latency (Negative Response)')

spikeGraphName = 'STALatencyIndivHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now we can also go unit by unit, and pull the first maximum that we
%encounter. WARNING WARNING STOP HERE ADJUST


%just plot averages


hFig = figure;
hold on
set(hFig, 'Position', [10 80 500 500])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
tester = mean(compTimingPos(:,dmrPV)');
plot([timeStep:timeStep:timeStep*100],tester(end:-1:1)/max(tester),'r')
tester = mean(compTimingPos(:,dmrMSN)');
plot([timeStep:timeStep:timeStep*100],tester(end:-1:1)/max(tester),'k')
tester = mean(compTimingNeg(:,dmrPV)');
plot([timeStep:timeStep:timeStep*100],-1*tester(end:-1:1)/min(tester),'r')
tester = mean(compTimingNeg(:,dmrMSN)');
plot([timeStep:timeStep:timeStep*100],-1*tester(end:-1:1)/min(tester),'k')

xlabel('Time (sec)')
ylabel('Normalized Stimulus Strength (a.u.)')
set(gca,'TickDir','out');

spikeGraphName = 'STALatencyPopulationMean';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% lets generate population average STAs
staFSI = zeros(size(reshape(bigSTAstore(1,:),length(faxis),[])));
for i = 1:length(dmrPV)
    tempStore = reshape(bigSTAstore(dmrPV(i),:),length(faxis),[]);
    tempStore = tempStore/max(max(bigSTAstore(dmrPV(i),:)));
    staFSI = staFSI + tempStore;
end

staMSN = zeros(size(reshape(bigSTAstore(1,:),length(faxis),[])));
for i = 1:length(dmrMSN)
    tempStore = reshape(bigSTAstore(dmrMSN(i),:),length(faxis),[]);
    tempStore = tempStore/max(max(bigSTAstore(dmrMSN(i),:)));
    staMSN = staMSN + tempStore;
end


staFSI = staFSI/max(max(staFSI));
staMSN = staMSN/max(max(staMSN));

figure
subplot(2,1,1)
clims = [0 1];
imagesc(staFSI,clims)
colorbar
set(gca, 'YDir', 'normal');
colormap(map)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('FSI Population STA')
subplot(2,1,2)
clims = [0 1];
imagesc(staMSN,clims)
colorbar
set(gca, 'YDir', 'normal');
colormap(map)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('MSN Population STA')


%% now lets pull out the RTFs
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

% for i = 1:length(dmrPV)
%     hFig = figure;
%     set(hFig, 'Position', [10 10 500 1000])
%     %make subplot
%     subplot(3,1,1)
%     imagesc(reshape(bigSTAstore(dmrPV(i),:),length(faxis),[]))
%     colormap('parula')
%     title('STA')
%     subplot(3,1,2)
%     imagesc(newSTA(:,:,dmrPV(i)))
%     colormap('parula')
%     title('STA Sig')
%     subplot(3,1,3)
%     imagesc(flipRTF(:,:,dmrPV(i)))
%     colormap('parula')
% end
% 
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


% hFig = figure;
% set(hFig, 'Position', [10 80 1900 1000])
% subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
% axisVal = ceil(sqrt(length(dmrPV)));
% for i = 1:length(dmrPV)
%     subplot(axisVal,axisVal,i)
%     imagesc(flipRTF(:,:,dmrPV(i)))
%     set(gca, 'YDir', 'normal');
%     colormap('parula')
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
% end
% 
% spikeGraphName = 'fsiRTFs';
% savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')
% 
% hFig = figure;
% set(hFig, 'Position', [10 80 1900 1000])
% subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
% axisVal = ceil(sqrt(length(dmrMSN)));
% for i = 1:length(dmrMSN)
%     subplot(axisVal,axisVal,i)
%     imagesc(flipRTF(:,:,dmrMSN(i)))
%     set(gca, 'YDir', 'normal');
%     colormap('parula')
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
% end
% 
% spikeGraphName = 'msnRTFs';
% savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')



%lets generate population average RTFs
rtfFSI = zeros(size(flipRTF(:,:,1)));
for i = 1:length(dmrPV)
    tempStore = flipRTF(:,:,dmrPV(i));
    tempStore = tempStore/max(max(tempStore));
    rtfFSI = rtfFSI + tempStore;
end

rtfMSN = zeros(size(flipRTF(:,:,1)));
for i = 1:length(dmrMSN)
    tempStore = flipRTF(:,:,dmrMSN(i));
    tempStore = tempStore/max(max(tempStore));
    rtfMSN = rtfMSN + tempStore;
end


rtfFSI = rtfFSI/max(max(rtfFSI));
rtfMSN = rtfMSN/max(max(rtfMSN));

figure
subplot(2,1,1)
imagesc(rtfFSI)
colorbar
set(gca, 'YDir', 'normal');
colormap(map)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('FSI Population RTF')
xlabel('TemporalMod')
ylabel('SpectralMod')
subplot(2,1,2)
imagesc(rtfMSN)
colorbar
set(gca, 'YDir', 'normal');
colormap(map)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('MSN Population RTF')
xlabel('TemporalMod')
ylabel('SpectralMod')


% %lets generate population average RTFs
% rtfFSI = zeros(size(rtfStore(:,:,1)));
% for i = 1:length(dmrPV)
%     tempStore = rtfStore(:,:,dmrPV(i));
%     tempStore = tempStore/max(max(tempStore));
%     rtfFSI = rtfFSI + tempStore;
% end
% 
% rtfMSN = zeros(size(rtfStore(:,:,1)));
% for i = 1:length(dmrMSN)
%     tempStore = rtfStore(:,:,dmrMSN(i));
%     tempStore = tempStore/max(max(tempStore));
%     rtfMSN = rtfMSN + tempStore;
% end
% 
% 
% rtfFSI = rtfFSI/max(max(rtfFSI));
% rtfMSN = rtfMSN/max(max(rtfMSN));
% 
% figure
% subplot(2,1,1)
% imagesc(rtfFSI)
% colorbar
% set(gca, 'YDir', 'normal');
% colormap('parula')
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% title('FSI Population RTF')
% subplot(2,1,2)
% imagesc(rtfMSN)
% colorbar
% set(gca, 'YDir', 'normal');
% colormap('parula')
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% title('MSN Population RTF')
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

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
% axisDef = ceil(sqrt(length(bigDBStore)));
% figure
% for i = 1:length(bigDBStore)
%     subplot(axisDef,axisDef,i)
%     if ismember(i,dmrPV)
%         plot(rtfModTemp(:,i),'r')
%     elseif ismember(i,dmrMSN)
%         plot(rtfModTemp(:,i),'k')
%     else
%         plot(rtfModTemp(:,i),'Color',[0.7 0.7 0.7])
%     end
% %     plot(rtfModTemp(:,i))
% end



% axisDef = ceil(sqrt(length(bigDBStore)));
% figure
% for i = 1:length(bigDBStore)
%     subplot(axisDef,axisDef,i)
%     if ismember(i,dmrPV)
%         plot(rtfModSpect(:,i),'r')
%     elseif ismember(i,dmrMSN)
%         plot(rtfModSpect(:,i),'k')
%     else
%         plot(rtfModSpect(:,i),'Color',[0.7 0.7 0.7])
%     end
% end

%plot out temporal
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
figure
plot(rtfModTemp(:,dmrPV)./max(rtfModTemp(:,dmrPV)))

figure
plot(rtfModTemp(:,dmrMSN)./max(rtfModTemp(:,dmrMSN)))

for i = 1:length(bigMaster)
    [peakTempModVal(i) peakTempModInd(i)] = max(rtfModTemp(:,i));
end

tmfFold = tmf(16:end);


%plot out spectral
figure
plot(rtfModSpect(:,dmrPV)./max(rtfModSpect(:,dmrPV)))

figure
hold on
for i = 1:length(dmrPV)
    plot(rtfModSpect(:,dmrPV(i))./max(rtfModSpect(:,dmrPV(i)))+i/4)
end

figure
plot(rtfModSpect(:,dmrMSN)./max(rtfModSpect(:,dmrMSN)))

for i = 1:length(bigMaster)
    [peakSpectModVal(i) peakSpectModInd(i)] = max(rtfModSpect(:,i));
end



hFig = figure;
set(hFig, 'Position', [10 80 1000 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.05 0.05]);
subplot(2,2,1)
hist(tmfFold(peakTempModInd(dmrPV)),[0:mean(diff(tmfFold)):tmfFold(end)])
set(gca,'TickDir','out')
xlabel('Temporal BMF (cyc/sec)')
ylabel('Number of Units')
title('FSI Temporal Best Mod Freq')
subplot(2,2,3)
hist(tmfFold(peakTempModInd(dmrMSN)),[0:mean(diff(tmfFold)):tmfFold(end)])
set(gca,'TickDir','out')
xlabel('Temporal BMF (cyc/sec)')
ylabel('Number of Units')
title('MSN Temporal Best Mod Freq')
subplot(2,2,2)
hist(xmf(peakSpectModInd(dmrPV)),[0:mean(diff(xmf)):xmf(end)])
set(gca,'TickDir','out')
xlabel('Spectral BMF (cyc/octave)')
ylabel('Number of Units')
title('FSI Spectral Best Mod Freq')
subplot(2,2,4)
hist(xmf(peakSpectModInd(dmrMSN)),[0:mean(diff(xmf)):xmf(end)])
set(gca,'TickDir','out')
xlabel('Spectral BMF (cyc/octave)')
ylabel('Number of Units')
title('MSN Spectral Best Mod Freq')


spikeGraphName = 'BMFs Temporal Spectral';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%% Lets now look at non-linearities
nl = bigpspk' .* bigpxspk; 
nl = nl ./ bigpx;

%non-normalized

figure
subplot(2,1,1)
hold on
plot(nl(dmrPV,:)')
plot(mean(nl(dmrPV,:)),'r','LineWidth',2)
subplot(2,1,2)
hold on
plot(nl(dmrMSN,:)')
plot(mean(nl(dmrMSN,:)),'r','LineWidth',2)

figure
subplot(2,1,1)
hold on
plot(nl(dmrPV,:)'./max(nl(dmrPV,:)'))
plot(mean(nl(dmrPV,:))/max(mean(nl(dmrPV,:))),'r','LineWidth',2)
subplot(2,1,2)
hold on
plot(nl(dmrMSN,:)'./max(nl(dmrMSN,:)'))
plot(mean(nl(dmrMSN,:))/max(mean(nl(dmrMSN,:))),'r','LineWidth',2)

%lets calculate skew?
skewVal = sum(((nl - mean(nl')').^3)./(size(nl,2)*std(nl').^3)',2);
hFig = figure;
set(hFig, 'Position', [10 80 500 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.05 0.05]);
subplot(2,1,1)
hist(skewVal(dmrPV),[0:0.2:3])
set(gca,'TickDir','out')
xlabel('Skew of Non-Linearity')
ylabel('Number of Units')
title('FSI Skew')
subplot(2,1,2)
hist(skewVal(dmrMSN),[0:0.2:3])
set(gca,'TickDir','out')
xlabel('Skew of Non-Linearity')
ylabel('Number of Units')
title('MSN Skew')


spikeGraphName = 'Skew Plot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Generate a summary figure

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.04 0.04], [0.04 0.04]);
set(hFig, 'Position', [80 80 1600 1200])

%Row 1, columns 1-2 plot STA
subplot(6,5,1)
imagesc(stimulus(:,15000:25000))
colormap('parula')
set(gca,'XTick',[])
set(gca,'YTick',[])



%Row 1, columns 3-6 examples!

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.01], [0.04 0.04], [0.04 0.04]);

exampleFSIs = [1,3,9,21];
exampleMSNs = [11,18,34,53];
for i = 1:4
    subplot(8,10,2+i)
    imagesc(reshape(bigSTASigstore(dmrMSN(exampleMSNs(i)),:),length(faxis),[]))
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
for i = 1:4
    subplot(8,10,12+i)
    imagesc(reshape(bigSTASigstore(dmrPV(exampleFSIs(i)),:),length(faxis),[]))
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.04 0.04], [0.04 0.04]);

%Row 1, column 7-10, rate, PLI, FSI, NonLineSkew

%Rate change
histVectSpikeNum = [-1:0.1:1];
rateModMSN = hist((dmrSpikeNum(dmrMSN)/(10*60)- bigMaster(dmrMSN,8)')./(dmrSpikeNum(dmrMSN)/(10*60)+ bigMaster(dmrMSN,8)'),histVectSpikeNum);
rateModFSI = hist((dmrSpikeNum(dmrPV)/(10*60)- bigMaster(dmrPV,8)')./(dmrSpikeNum(dmrPV)/(10*60)+ bigMaster(dmrPV,8)'),histVectSpikeNum);
ylimRateMod = max([max(rateModMSN),max(rateModFSI)]);

subplot(8,10,7)
bar(histVectSpikeNum,rateModMSN,'k')
ylim([0 ylimRateMod])
xlim([-1.05 1.05])
set(gca,'TickDir','out')
set(gca,'xticklabel',[])
ylabel('#units')

subplot(8,10,17)
bar(histVectSpikeNum,rateModFSI,'r')
ylim([0 ylimRateMod])
xlim([-1.05 1.05])
set(gca,'TickDir','out')
ylabel('#units')
xlabel('Modulation Index')

%Phase Locking Index
histVectSpikeNum = [0:0.05:1];
pliMSN = hist(PLI(dmrMSN),histVectSpikeNum);
pliFSI = hist(PLI(dmrPV),histVectSpikeNum);
ylimPLI  = max([max(pliMSN),max(pliFSI)]);

subplot(8,10,8)
bar(histVectSpikeNum,pliMSN,'k')
ylim([0 ylimPLI])
xlim([histVectSpikeNum(1) histVectSpikeNum(end)])
set(gca,'TickDir','out')
set(gca,'xticklabel',[])
ylabel('#units')

subplot(8,10,18)
bar(histVectSpikeNum,pliFSI,'r')
ylim([0 ylimPLI])
xlim([histVectSpikeNum(1) histVectSpikeNum(end)])
set(gca,'TickDir','out')
ylabel('#units')
xlabel('Phase Lock Index')

%FSI (feature selectivity index)


%Skew (nonlinearity)
skewVect = [0:0.2:3];
skewMSN = hist(skewVal(dmrMSN),skewVect);
skewFSI = hist(skewVal(dmrPV),skewVect);
ylimSkew = max([max(skewMSN),max(skewFSI)]);

subplot(8,10,10)
bar(skewVect,skewMSN,'k')
set(gca,'TickDir','out')
set(gca,'xticklabel',[])
ylim([0 ylimSkew])
xlim([skewVect(1) skewVect(end)])
% xlabel('Skew of Non-Linearity')
ylabel('Number of Units')
title('MSN Skew')

subplot(8,10,20)
bar(skewVect,skewFSI,'k')
set(gca,'TickDir','out')
% set(gca,'xticklabel',[])
ylim([0 ylimSkew])
xlim([skewVect(1) skewVect(end)])
% xlabel('Skew of Non-Linearity')
ylabel('Number of Units')
title('FSI Skew')



%Row 2, columns 1, 2 plot tuning curve and STA

subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.04 0.04], [0.04 0.04]);
%last FSI is a pretty good example
subplot(4,6,7)
testData = binValBigStore(:,:,dmrMSN(55));
%we need to stack in spaces for frequencies not represented
testData(17:21,:) = 0;
imagesc(testData)
colormap('parula')
set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(4,6,8)
imagesc(newSTA(:,:,dmrMSN(55)))
colormap('parula')
set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(4,6,13)
testData = binValBigStore(:,:,328);
%we need to stack in spaces for frequencies not represented
testData(17:21,:) = 0;
imagesc(testData)
colormap('parula')
set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(4,6,14)
imagesc(newSTA(:,:,328))
colormap('parula')
set(gca,'XTick',[])
set(gca,'YTick',[])


subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.04 0.04], [0.04 0.04]);

%Row 2 column 3:  plot pie charts of significant tone response vs STA response. 
indPosSig = 9;
indNegSig = 10;
holder = bigMaster(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
detMSN = holder(:,1) + holder(:,2);
detMSNHist = hist(detMSN,[-2:1:1]);

% pie(detMSNHist)
% labels = {'Neg','Mix','None','Pos'};
% detZero = find(det == 0);
% labels(detZero) = [];

holder = bigMaster(findPVs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
detFSI = holder(:,1) + holder(:,2);
detFSIHist = hist(detFSI,[-2:1:1]);

% pie(detFSIHist)
% labels = {'Neg','Mix','None','Pos'};
% detZero = find(det == 0);
% labels(detZero) = [];

%basically, just going to compare for non-zero values. 
basicDetMSN = detMSN;
basicDetMSN(basicDetMSN ~=0) = 1;

basicDetFSI = detFSI;
basicDetFSI(basicDetFSI ~=0) = 1;


justDMRMSN = findMSNs;
justDMRMSN(ismember(justDMRMSN,fewSpikes)) = [];
justDMRMSN = intersect(justDMRMSN,sigUnitPrct);
[C ia justDMRMSN] = intersect(justDMRMSN,findMSNs);

justDMRFSI= findPVs;
justDMRFSI(ismember(justDMRFSI,fewSpikes)) = [];
justDMRFSI = intersect(justDMRFSI,sigUnitPrct);
[C ia justDMRFSI] = intersect(justDMRFSI,findPVs);

combDetMSN = basicDetMSN;
combDetMSN(justDMRMSN) = combDetMSN(justDMRMSN) + 2;
combDetHistMSN = hist(combDetMSN,[0:3]);

combDetFSI = basicDetFSI;
combDetFSI(justDMRFSI) = combDetFSI(justDMRFSI) + 2;
combDetHistFSI = hist(combDetFSI,[0:3]);

subplot(4,6,9)
pie(combDetHistMSN)

subplot(4,6,15)
pie(combDetHistFSI)

subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.03], [0.04 0.04], [0.04 0.04]);
%Row 2, column 4-6, now we make plots for width/latency/bf?
subplot(4,6,10)
hold on
maxVal = max(max(widthCompMSN));
plot(widthCompMSN(2,:),widthCompMSN(1,:),'k.')
plot([0 maxVal],[0 maxVal],'Color',[0.7 0.7 0.7])
set(gca,'TickDir','out')
set(gca,'xticklabel',[])
xlim([0 maxVal])
ylim([0 maxVal])
axis square


subplot(4,6,16)
hold on
maxVal = max(max(widthCompFSI));
plot(widthCompFSI(2,:),widthCompFSI(1,:),'r.')
plot([0 maxVal],[0 maxVal],'Color',[0.7 0.7 0.7])
set(gca,'TickDir','out')
xlim([0 maxVal])
ylim([0 maxVal])
axis square


%plot BFs
maxVal = 32000;
subplot(4,6,11)
hold on
loglog(staBFPos(dmrMSN),toneBF(dmrMSN),'k.')
plot([4000 64000],[4000 64000],'Color',[0.7 0.7 0.7])
set(gca,'TickDir','out')
set(gca,'xticklabel',[])
% xlabel('STA BF (Hz)')
% ylabel('Pure Tone BF (Hz)')
xlim([0 maxVal])
ylim([0 maxVal])
axis square


subplot(4,6,17)
hold on
loglog(staBFPos(dmrPV),toneBF(dmrPV),'r.')
plot([4000 64000],[4000 64000],'Color',[0.7 0.7 0.7])
set(gca,'TickDir','out')
xlabel('STA BF (Hz)')
ylabel('Pure Tone BF (Hz)')
xlim([0 maxVal])
ylim([0 maxVal])
axis square

%plot latencies

subplot(4,6,12)
hold on
maxVal = max([max(peakTimePos(dmrMSN)) max(latStore(dmrMSN))]);
plot(peakTimePos(dmrMSN),latStore(dmrMSN),'k.')
plot([0 maxVal],[0 maxVal],'Color',[0.7 0.7 0.7])
set(gca,'TickDir','out')
ylabel('ToneLat')
xlim([0 maxVal])
ylim([0 maxVal])
axis square

subplot(4,6,18)
hold on
maxVal = max([max(peakTimePos(dmrPV)) max(latStore(dmrPV))]);
plot(peakTimePos(dmrPV),latStore(dmrPV),'r.')
plot([0 maxVal],[0 maxVal],'Color',[0.7 0.7 0.7])
set(gca,'TickDir','out')
ylabel('ToneLat')
xlim([0 maxVal])
ylim([0 maxVal])
axis square


%Row 3, columns 1 and 2: STA to RTF transformation example. Again, use last
%dmrPV
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.04 0.04], [0.04 0.04]);
subplot(4,8,25)
imagesc(newSTA(:,:,328))
colormap('parula')
set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(4,8,26)
imagesc(flipRTF(:,:,dmrPV(end)))
set(gca, 'YDir', 'normal');
colormap('parula')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'TickDir','out')

%Row 3 columns 3 and 4: example RTFs
targetMSNRTF = [5,42];
targetFSIRTF = [16,18];
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.04 0.04], [0.04 0.04]);
for i = 1:2
    subplot(8,8,50+i)
    imagesc(flipRTF(:,:,dmrMSN(targetMSNRTF(i))))
    set(gca, 'YDir', 'normal');
    colormap('parula')
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'TickDir','out')
    subplot(8,8,58+i)
    imagesc(flipRTF(:,:,dmrPV(targetFSIRTF(i))))
    set(gca, 'YDir', 'normal');
    colormap('parula')
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'TickDir','out')
    
end

%Row 3 column 5-6, plot out population RTFs

subplot(4,8,29)
imagesc(rtfMSN)
colorbar
set(gca, 'YDir', 'normal');
colormap('parula')
% set(gca,'xtick',[])
set(gca,'xticklabel',[])
% set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('MSN Population RTF')

subplot(4,8,30)
imagesc(rtfFSI)
colorbar
set(gca, 'YDir', 'normal');
colormap('parula')
% set(gca,'xtick',[])
set(gca,'xticklabel',[])
% set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('FSI Population RTF')

%Row 3, column 7-8

tempModHistMSN = hist(tmfFold(peakTempModInd(dmrMSN)),[0:mean(diff(tmfFold)):tmfFold(end)]);
tempModHistFSI = hist(tmfFold(peakTempModInd(dmrPV)),[0:mean(diff(tmfFold)):tmfFold(end)]);
ylimTempMod = max([max(tempModHistMSN),max(tempModHistFSI)]);

subplot(8,8,55)
bar([0:mean(diff(tmfFold)):tmfFold(end)],tempModHistMSN,'k')
set(gca,'TickDir','out')
% xlabel('Temporal BMF (cyc/sec)')
ylabel('Number of Units')
ylim([0 ylimTempMod])
title('MSN Temporal Best Mod Freq')

subplot(8,8,63)
bar([0:mean(diff(tmfFold)):tmfFold(end)],tempModHistFSI,'r')
set(gca,'TickDir','out')
xlabel('Temporal BMF (cyc/sec)')
ylabel('Number of Units')
ylim([0 ylimTempMod])
title('FSI Temporal Best Mod Freq')

spectModHistMSN = hist(xmf(peakSpectModInd(dmrMSN)),[0:mean(diff(xmf)):xmf(end)]);
spectModHistFSI = hist(xmf(peakSpectModInd(dmrPV)),[0:mean(diff(xmf)):xmf(end)]);
ylimSpectMod = max([max(spectModHistMSN),max(spectModHistFSI)]);

subplot(8,8,56)
bar([0:mean(diff(xmf)):xmf(end)],spectModHistMSN,'k')
set(gca,'TickDir','out')
ylim([0 ylimSpectMod])
% xlabel('Spectral BMF (cyc/octave)')
% ylabel('Number of Units')
title('MSN Spectral Best Mod Freq')

subplot(8,8,64)
bar([0:mean(diff(xmf)):xmf(end)],spectModHistFSI,'r')
set(gca,'TickDir','out')
ylim([0 ylimSpectMod])
xlabel('Spectral BMF (cyc/octave)')
% ylabel('Number of Units')
title('FSI Spectral Best Mod Freq')



spikeGraphName = 'DMRSummaryPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
print(hFig,spikeGraphName,'-deps','-r0')














