
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

%plot out spiking rate for DMR and PLI

%calculate phase locking index
% PLI = (max(bigSTAstore') - min(bigSTAstore'))./(dmrSpikeNum/(10*60)*sqrt(8));
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


%plot out STAs

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
axisVal = ceil(sqrt(length(dmrPV)));
for i = 1:length(dmrPV)
    subplot(axisVal,axisVal,i)
    imagesc(reshape(bigSTASigstore(dmrPV(i),:),length(faxis),[]))
    colormap('parula')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
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
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
axisVal = ceil(sqrt(length(dmrMSN)));
for i = 1:length(dmrMSN)
    subplot(axisVal,axisVal,i)
    imagesc(reshape(bigSTASigstore(dmrMSN(i),:),length(faxis),[]))
    colormap('parula')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
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

%to determine peak and width, lets compress along time axis. 
smoothWind = 7;
% smoothWind = 1;
for i = 1:length(bigDBStore)
    compTimingPos(:,i) = sum(posSTA(:,:,i));
    compTimingNeg(:,i) =  sum(negSTA(:,:,i));
    compWidthPos(:,i) = smooth(sum(posSTA(:,:,i)'),smoothWind);
    compWidthNeg(:,i) = smooth(sum(negSTA(:,:,i)'),smoothWind);
    
end



%Lets examine BFs
[maxVal staBFPos] = max(compWidthPos);
staBFPos = faxis(staBFPos);

[maxVal staBFNeg] = min(compWidthNeg);
staBFNeg = faxis(staBFNeg);

%BF store from master data is sketchy. Pull directly from final amplitude
%of binDiff. 

[maxVal toneBF] = max(squeeze(binValBigStore(:,end,:)));
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
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
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
    plot(toneFreqs,toneWidth/max(toneWidth),'k')
    plot(faxis,compWidthPos(:,dmrMSN(i))/max(compWidthPos(:,dmrMSN(i))),'r')
    plot(faxis,-1*compWidthNeg(:,dmrMSN(i))/min(compWidthNeg(:,dmrMSN(i))),'b')
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


%Now lets go to 0.2 of height
%now lets calculate widths based on half height. 
heightBench = 0.2;
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
title('Change in FSI 20% Width')
subplot(2,1,2)
hold on
for i = 1:length(dmrMSN)
    plot([1 2],[widthCompMSN(1,i) widthCompMSN(2,i)],'Color',[0.7 0.7 0.7])
end
errorbar([1 2],[nanmean(widthCompMSN(1,:)) nanmean(widthCompMSN(2,:))],[nanstd(widthCompMSN(1,:))/sqrt(length(dmrMSN)) nanstd(widthCompMSN(2,:))/sqrt(length(dmrMSN))],'k','LineWidth',2)
xlim([0.5 2.5])
set(gca,'TickDir','out');
title('Change in MSN 20% Width')

spikeGraphName = 'Width20%STAvsPureTone';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')






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


hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
axisVal = ceil(sqrt(length(dmrPV)));
for i = 1:length(dmrPV)
    subplot(axisVal,axisVal,i)
    imagesc(flipRTF(:,:,dmrPV(i)))
    set(gca, 'YDir', 'normal');
    colormap('parula')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

spikeGraphName = 'fsiRTFs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.05 0.05]);
axisVal = ceil(sqrt(length(dmrMSN)));
for i = 1:length(dmrMSN)
    subplot(axisVal,axisVal,i)
    imagesc(flipRTF(:,:,dmrMSN(i)))
    set(gca, 'YDir', 'normal');
    colormap('parula')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

spikeGraphName = 'msnRTFs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

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



