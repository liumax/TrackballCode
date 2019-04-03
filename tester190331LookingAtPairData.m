%This code is janky code to look at identified pairs. The goal is the
%produce plots of the cross correlograms (without smoothing),
%representations of responses, and things like difference in best
%frequency.



targetFiles = {'180315_ML180306C_R17_3218mid1_fullTuningAnalysis.mat',...
    '180718_ML180619B_L_AudStr_pen1_3000_fullTuningAnalysis.mat',...
    '180718_ML180619C_R_AudStr_pen1_3000_fullTuningAnalysis.mat',...
    '190123ML181105E_RAudStr3526pen1rec1tuningAndDMRAnalysis.mat',...
    '190205ML181105C_RAudStr3667pen2rec1tuningDMRAnalysis.mat',...
    '190206ML181105F_RAudStr3633pen2rec1tuningDMRAnalysis.mat'};


orderFile = [1	106	6
2	257	166
2	257	27
2	257	32
2	269	265
2	270	32
2	270	199
3	25	256
3	25	112
3	25	258
4	41	70
4	41	150
4	271	278
4	274	278
5	284	283
5	292	7
6	265	270
6	273	272
6	273	277
6	273	84
6	273	40
6	283	15
6	283	40
6	283	95];

overCount = 1;
truexCorr = zeros(41,length(orderFile));
widexCorrFSI = zeros(41,length(orderFile));
widexCorrMSN = zeros(41,length(orderFile));
rangeBounds = zeros(41,length(orderFile));
rangeBoundsPRCTL = zeros(41,2,length(orderFile));
toneStore = zeros(41,length(orderFile));
nonToneStore = zeros(41,length(orderFile));
remStore = zeros(41,length(orderFile));

numShuff = 1000;

for i = 1:length(targetFiles)
    disp(targetFiles{i})
    load(targetFiles{i});
    %first, figure out how many pairs I'm looking for
    tarUnits = find(orderFile(:,1) == i);
    numUnits = length(tarUnits);
    
    tester = s.DesignationName;
    numName = zeros(length(tester),1);
    for j = 1:length(tester);
        numName(j) = str2num(tester{j}(4:end));
    end
    fsiInd = [];
    msnInd = [];
    shankInd = [];
    %determine number of tone reps and number of DBs, then subtract
        %these value from the raster.
    numTones = s.SoundData.ToneRepetitions;
    numDBs = s.SoundData.NumDBs;
    subNum = numTones*numDBs*s.SoundData.WhiteNoise
    
    for j = 1:numUnits
        %first find FSI
        fsiInd(j) = find(numName == orderFile(tarUnits(j),2));
        %then find MSN
        msnInd(j) = find(numName == orderFile(tarUnits(j),3));
        %determine BF store
        [bfInd] = functionCellStringFind(masterHeader,'BF');
        fsiBFStore(overCount + j - 1) = masterData(fsiInd(j),bfInd);
        msnBFStore(overCount + j - 1) = masterData(msnInd(j),bfInd);
        %store the frequency organized histograms, maybe this will be a
        %nice way to visualize interaction?
        fsiHistStore{overCount + j - 1} = s.(tester{fsiInd(j)}).FrequencyHistograms(1+s.SoundData.WhiteNoise:end,:);
        msnHistStore{overCount + j - 1} = s.(tester{msnInd(j)}).FrequencyHistograms(1+s.SoundData.WhiteNoise:end,:);
        histBinVector{overCount + j - 1} = s.(tester{fsiInd(j)}).HistBinVector;
        %then find shank value based on peakChan
        shankInd(j) = s.PeakChanVals(fsiInd(j));
        %also, lets grab rasters and shit
        
        tmpRasters = s.(tester{fsiInd(j)}).AllRasters;
        tmpRasters(tmpRasters(:,3) <= subNum,:) = [];
        rastersFSI{overCount + j - 1} = tmpRasters;
        tmpRasters = s.(tester{msnInd(j)}).AllRasters;
        tmpRasters(tmpRasters(:,3) <= subNum,:) = [];
        rastersMSN{overCount + j - 1} = tmpRasters;
        soundDataStore(overCount + j - 1,1) = length(s.SoundData.UniqueFrequencies);
        soundFreqs{overCount + j - 1} = s.SoundData.UniqueFrequencies;
        %grab the tone period data so i can pull BFs and width etc. 
        toneWind = 2;
        binDiffValFSI{overCount + j - 1} = squeeze(s.(tester{fsiInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind));
        binDiffValMSN{overCount + j - 1} = squeeze(s.(tester{msnInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind));
        
        
        %and grab correlation coefficient of binned baseline subtracted tone window responses
        toneWind = 2;
        val1 = reshape(squeeze(s.(tester{fsiInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind)),1,[]);
        val2 = reshape(squeeze(s.(tester{msnInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind)),1,[]);
        tmpCorr = corrcoef(val1,val2);
        binRespCorrTone(overCount + j - 1) = tmpCorr(2);
        %now calculate significance
        for k = 1:numShuff
            shuffVal = randperm(length(val2));
            shuffVal = val2(shuffVal);
            tempCorrVal = corrcoef(val1,shuffVal);
            shuffValStoreTone(k,overCount + j - 1) = tempCorrVal(2);
        end
        %also get onset. 
        toneWind = 1;
        val1 = reshape(squeeze(s.(tester{fsiInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind)),1,[]);
        val2 = reshape(squeeze(s.(tester{msnInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind)),1,[]);
        tmpCorr = corrcoef(val1,val2);
        binRespCorrFast(overCount + j - 1) = tmpCorr(2);
        %now calculate significance
        for k = 1:numShuff
            shuffVal = randperm(length(val2));
            shuffVal = val2(shuffVal);
            tempCorrVal = corrcoef(val1,shuffVal);
            shuffValStoreFast(k,overCount + j - 1) = tempCorrVal(2);
        end
        %also get Long. 
        toneWind = 3;
        val1 = reshape(squeeze(s.(tester{fsiInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind)),1,[]);
        val2 = reshape(squeeze(s.(tester{msnInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,toneWind)),1,[]);
        tmpCorr = corrcoef(val1,val2);
        binRespCorrSlow(overCount + j - 1) = tmpCorr(2);
        %now calculate significance
        for k = 1:numShuff
            shuffVal = randperm(length(val2));
            shuffVal = val2(shuffVal);
            tempCorrVal = corrcoef(val1,shuffVal);
            shuffValStoreSlow(k,overCount + j - 1) = tempCorrVal(2);
        end
    end
    %correct shank ind to report shank number
    shankInd(shankInd <= 32) = 1;
    shankInd(shankInd > 1) = 2;
    %also need to determine new indexes for each shank, since organized by
    %shank for correlation data. 
    if i == 1
        shank1order = sort(numName(s.PeakChanVals <= 32));
    else
        shank1order = sort(numName(s.PeakChanVals < 32));
    end
    shank2order = sort(numName(s.PeakChanVals > 32));
    
    
    %NOW lets go and look for the cross correlogram. 
    
    for j = 1:numUnits
        trueName = strcat('trueStoreShank',num2str(shankInd(j)));
        wideName = strcat('wideStoreShank',num2str(shankInd(j)));
        shuffName = strcat('ShuffStore',num2str(shankInd(j)));
        toneName = strcat('Tone',num2str(shankInd(j)));
        nonToneName = strcat('NonTone',num2str(shankInd(j)));
        remName = strcat('RemStore',num2str(shankInd(j)));
        if shankInd(j) == 1
            newFSI = find(shank1order == numName(fsiInd(j)));
            newMSN = find(shank1order == numName(msnInd(j)));
        elseif shankInd(j) == 2
            newFSI = find(shank2order == numName(fsiInd(j)));
            newMSN = find(shank2order == numName(msnInd(j)));
        end
        truexCorr(:,overCount + j - 1) = squeeze(corrData.(trueName)(newFSI,newMSN,:));
        widexCorrFSI(:,overCount + j - 1) = squeeze(corrData.(wideName)(newFSI,newFSI,:));
        widexCorrMSN(:,overCount + j - 1) = squeeze(corrData.(wideName)(newMSN,newMSN,:));
        toneStore(:,overCount + j - 1) = squeeze(corrData.(toneName)(newFSI,newMSN,:));
        nonToneStore(:,overCount + j - 1) = squeeze(corrData.(nonToneName)(newFSI,newMSN,:));
        remStore(:,overCount + j - 1) = squeeze(corrData.(remName)(newFSI,newMSN,:));
        rangeBounds(:,overCount + j - 1) = mean(corrData.(shuffName){newFSI,newMSN}');
        rangeBoundsPRCTL(:,:,overCount + j - 1) = prctile(corrData.(shuffName){newFSI,newMSN},[0.5 99.5],2);
        
        recNum(overCount + j - 1) = i;
    end
    
    overCount = overCount + numUnits;
    
end


%pull some other characteristics
overCount = 1;
sigCut = 0.01;
for i = 1:length(targetFiles)
    disp(targetFiles{i})
    load(targetFiles{i});
    %first, figure out how many pairs I'm looking for
    tarUnits = find(orderFile(:,1) == i);
    numUnits = length(tarUnits);
    
    tester = s.DesignationName;
    numName = zeros(length(tester),1);
    for j = 1:length(tester);
        numName(j) = str2num(tester{j}(4:end));
    end
    for j = 1:numUnits
        %first find FSI
        fsiInd(j) = find(numName == orderFile(tarUnits(j),2));
        %then find MSN
        msnInd(j) = find(numName == orderFile(tarUnits(j),3));
        
        %determine number of significant positive responses.
        
        %FSIs
        binResp = squeeze(s.(tester{fsiInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,2));
        sigResp = squeeze(s.(tester{fsiInd(j)}).BinSigVals(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,2));
        signResp = sign(binResp);
        signResp(signResp < 1) = NaN;
        signSigResp = sigResp.*signResp;
        signStore{overCount + j -1,1} = signSigResp;
        findSig = length(find(signSigResp < sigCut));
        sigNum(overCount + j - 1,1) = findSig;
        %now lets pull the responses from binDiffTone. 
%         cleanBinFSI{overCount + j - 1} = zeros(size(binDiffValFSI{overCount + j -1}));
        cleanBinFSI{overCount + j - 1} = binDiffValFSI{overCount + j -1}.*signResp;
        %MSNs
        binResp = squeeze(s.(tester{msnInd(j)}).BinDiff(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,2));
        sigResp = squeeze(s.(tester{msnInd(j)}).BinSigVals(1+s.SoundData.WhiteNoise:end,1:s.SoundData.NumDBs,2));
        signResp = sign(binResp);
        signResp(signResp < 1) = NaN;
        signSigResp = sigResp.*signResp;
        signStore{overCount + j -1,2} = signSigResp;
        findSig = length(find(signSigResp < sigCut & signSigResp > 0));
        sigNum(overCount + j - 1,2) = findSig;
        %now lets pull the responses from binDiffTone. 
%         cleanBinMSN{overCount + j - 1} = zeros(size(binDiffValMSN{overCount + j -1}));
        cleanBinMSN{overCount + j - 1} = binDiffValMSN{overCount + j -1}.*signResp;
        
        %also store latency map.
        latMapFSI{overCount + j - 1} = s.(tester{fsiInd(j)}).LatencyMap(1+s.SoundData.WhiteNoise:end,:);
        latMapMSN{overCount + j - 1} = s.(tester{msnInd(j)}).LatencyMap(1+s.SoundData.WhiteNoise:end,:);
        
        
        
    end
    overCount = overCount + j;
end

%now lets find things!
counter = 1;
findBad = [];
minVal = 5;
for i = 1:2
    tmpVals = sigNum(:,i);
    tmpVals = find(tmpVals < minVal);
    findBad(counter:counter+length(tmpVals)-1) = tmpVals;
    counter = counter + length(tmpVals);
end

findBad = unique(findBad);

%make good index
findGood = [1:length(orderFile)];
findGood(findBad) = [];
    

%plot out true correlations
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
    subplot(6,4,i)
    hold on
    plot(xAxVals,truexCorr(:,i),'r','LineWidth',2)
    plot(xAxVals,rangeBounds(:,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,2,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,1,i),'Color',[0.7 0.7 0.7])
end


spikeGraphName = 'AllCrossCorrSig';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%plot out tone vs non-tone
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
    subplot(6,4,i)
    hold on
    plot(xAxVals,toneStore(:,i),'m','LineWidth',2)
    plot(xAxVals,nonToneStore(:,i),'k','LineWidth',2)
end


spikeGraphName = 'AllCrossToneVsNonTone';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%plot out correlations with synchronous spikes removed
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
    subplot(6,4,i)
    hold on
    plot(xAxVals,remStore(:,i),'k','LineWidth',2)
    plot(xAxVals,rangeBounds(:,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,2,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,1,i),'Color',[0.7 0.7 0.7])
end


spikeGraphName = 'AllCrossCorrSigSynchRemoved';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% lets look at tuning curve correlations
figure
hist(binRespCorrTone,[-1:0.1:1])

figure
hist(binRespCorrFast,[-1:0.1:1])

figure
hist(binRespCorrSlow,[-1:0.1:1])


% figure
% hold on
% for i = 1:length(orderFile)
%     plot([binRespCorrTone(i),binRespCorrFast(i),binRespCorrSlow(i)])
% end

%lets calculate which are significantly + or -
findCorrSig = zeros(length(orderFile),1);
for i = 1:length(orderFile)
    bounds = prctile(shuffValStoreTone(:,i),[0.5 99.5]);
    if binRespCorrTone(i) < bounds(1) || binRespCorrTone(i) > bounds(2)
        findCorrSig(i) = 1;
    end
end

histVectCorr = [-1:0.1:1];
sketchVect = [4,6,9,11,12,13];
allBinToneCorr = hist(binRespCorrTone,histVectCorr);
sigBinToneCorr = hist(binRespCorrTone(findCorrSig==1),histVectCorr);
selBinToneCorr = hist(binRespCorrTone(intersect(find(findCorrSig == 1),sketchVect)),histVectCorr);
justSketch = hist(binRespCorrTone(sketchVect),histVectCorr);

figure
hold on
bar(histVectCorr,allBinToneCorr,'w')
bar(histVectCorr,sigBinToneCorr,'k')
bar(histVectCorr,justSketch,'c')
bar(histVectCorr,selBinToneCorr,'r')

% 
% %lets now plot out the 70dB response curves for FSIs and MSNs
% for j = 1:5
%     hFig = figure;
%     set(hFig, 'Position', [10 80 1800 1200])
%     borderControl = ceil(sqrt(length(orderFile)));
%     subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
%     % xAxVals = [-0.01:0.0005:0.01];
%     for i = 1:length(orderFile);
%         subplot(borderControl,borderControl,i)
%     %     subplot(6,4,i)
%         hold on
%         plot(binDiffValFSI{i}(:,end-j+1),'r','LineWidth',2)
%         plot(binDiffValMSN{i}(:,end-j+1),'k','LineWidth',2)
%     end
% end
% 
%now lets plot out the histograms for each. 
for i = 1:length(orderFile)
    hFig = figure;
    set(hFig, 'Position', [10 80 500 1000])
%     borderControl = ceil(sqrt(length(orderFile)));
    subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
    numFreqs = size(fsiHistStore{i},1);
    for j = 1:numFreqs
        hold on
        plot(histBinVector{i},fsiHistStore{(i)}(j,:)/max(max(fsiHistStore{(i)}))+j,'r')
        plot(histBinVector{i},msnHistStore{(i)}(j,:)/max(max(msnHistStore{(i)}))+j,'k')
        plot(histBinVector{i},zeros(length(fsiHistStore{(i)}(j,:)),1)+j,'g')
    end
    plot([0 0],[1 numFreqs + 1],'Color',[0.7 0.7 0.7],'LineWidth',2)
    xlim([-0.1 0.2])
    ylim([1 numFreqs+1])
    set(gca,'TickDir','out')
    spikeGraphName = strcat('FrequencyHist-',num2str(i));
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end

%% Lets plot Rasters!

for i = 1:length(orderFile)
    hFig = figure;
    set(hFig, 'Position', [10 80 1000 1000])
    subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
    subplot(1,2,1)
    plot(rastersFSI{i}(:,1),rastersFSI{i}(:,3),'k.')
    subplot(1,2,2)
    plot(rastersMSN{i}(:,1),rastersMSN{i}(:,3),'k.')
    spikeGraphName = strcat('RastersSet-',num2str(i));
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

%% Now lets try and pull more information about responses! 


%Widths!

%calculate widths with height based system
bigWidthHeightStore = [];
widthFrac = 0.5;
for i = 1:length(orderFile)
    testValues = squeeze(binDiffValFSI{i});
    for j = 1:size(testValues,2)
        [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues(:,j),widthFrac,3);
        bigWidthHeightStore(j,:,i) = widthVals;
        bigWidthMaxPosStore(j,i) = maxPos;
        bigWidthMaxValStore(j,i) = maxVal;
        bigWidthCutVal(j,i) = cutVal;
    end
end

%now lets just pull out the outward search, since this looks to be far more
%accurate
bigWidthSelWidth = bigWidthHeightStore(:,[1,3],:);

%MSNs
bigWidthHeightStoreMSN = [];
for i = 1:length(orderFile)
    testValues = squeeze(binDiffValMSN{i});
    for j = 1:size(testValues,2)
        [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(testValues(:,j),widthFrac,3);
        bigWidthHeightStoreMSN(j,:,i) = widthVals;
        bigWidthMaxPosStoreMSN(j,i) = maxPos;
        bigWidthMaxValStoreMSN(j,i) = maxVal;
        bigWidthCutValMSN(j,i) = cutVal;
    end
end

%now lets just pull out the outward search, since this looks to be far more
%accurate
bigWidthSelWidthMSN = bigWidthHeightStoreMSN(:,[1,3],:);

%now we generate the full width using both sides
tarCellHeightWidth = bigWidthSelWidth;

% pvMaxPos = bigWidthMaxPosStore(:,tarCells(tarPVs));
% msnMaxPos = bigWidthMaxPosStore(:,tarCells(tarMSNs));
tarCellMaxPos = bigWidthMaxPosStore;

%now lets determine whether significant values exist at each DB level.
sigCut = 0.01;
% abbrevSig = signSigStore(:,:,tarCells);
for i = 1:length(orderFile)
    tester = min(signStore{i,1});
    tester(tester <= sigCut) = 10;
    tester(tester < 10) = 0;
    tester(tester == 10) = 1;
    tester(isnan(tester)) = 0;
    try
        threshAmp(i,1) = find(tester == 1,1,'first');
        threshAmp(i,2) = length(tester);
    catch
        threshAmp(i,1) = NaN;
        threshAmp(i,2) = length(tester);
    end
    sigStoreBinary{i} = tester;
end

%ok, now we have significant values, max points, and the detected width in
%each direction. Lets try and pull things out now?

sigWidthStore = NaN(6,length(orderFile));
sideWarnLow = NaN(6,length(orderFile));
sideWarnHi = NaN(6,length(orderFile));
for i = 1:length(orderFile)
    for j = 1:length(sigStoreBinary{i})
        tmpData = sigStoreBinary{i};
        %first, check if there is a significant value
        if tmpData(j) == 1
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
            sigWidthStore(j+(6-threshAmp(i,2)),i) = hiWidth + lowWidth;
        end
    end
end


%now lets just pull out the outward search, since this looks to be far more
%accurate
bigWidthSelWidthMSN = bigWidthHeightStoreMSN(:,[1,3],:);

%now we generate the full width using both sides
tarCellHeightWidth = bigWidthSelWidthMSN;

% pvMaxPos = bigWidthMaxPosStore(:,tarCells(tarPVs));
% msnMaxPos = bigWidthMaxPosStore(:,tarCells(tarMSNs));
tarCellMaxPos = bigWidthMaxPosStoreMSN;

%now lets determine whether significant values exist at each DB level.
sigCut = 0.01;
% abbrevSig = signSigStore(:,:,tarCells);
for i = 1:length(orderFile)
    tester = min(signStore{i,2});
    tester(tester <= sigCut) = 10;
    tester(tester < 10) = 0;
    tester(tester == 10) = 1;
    tester(isnan(tester)) = 0;
    try
        threshAmpMSN(i,1) = find(tester == 1,1,'first');
        threshAmpMSN(i,2) = length(tester);
    catch
        threshAmpMSN(i,1) = NaN;
        threshAmpMSN(i,2) = length(tester);
    end
    sigStoreBinaryMSN{i} = tester;
end

%ok, now we have significant values, max points, and the detected width in
%each direction. Lets try and pull things out now?

sigWidthStoreMSN = NaN(6,length(orderFile));
sideWarnLowMSN = NaN(6,length(orderFile));
sideWarnHiMSN = NaN(6,length(orderFile));
for i = 1:length(orderFile)
    for j = 1:length(sigStoreBinaryMSN{i})
        tmpData = sigStoreBinaryMSN{i};
        %first, check if there is a significant value
        if tmpData(j) == 1
            %now check for width values. If NaN on either side, determine
            %edge.
            %start with lower edge. 
            if isnan(tarCellHeightWidth(j,1,i))
                lowWidth = tarCellMaxPos(j,i);
                sideWarnLowMSN(j,i) = 1;
            else
                lowWidth = tarCellHeightWidth(j,1,i);
            end
            if isnan(tarCellHeightWidth(j,2,i))
                hiWidth = 16 - tarCellMaxPos(j,i);
                sideWarnHiMSN(j,i) = 1;
            else
                hiWidth = tarCellHeightWidth(j,2,i);
            end
            sigWidthStoreMSN(j+(6-threshAmpMSN(i,2)),i) = hiWidth + lowWidth;
        end
    end
end

%remove bad ones
sigWidthStoreSel = sigWidthStore;
sigWidthStoreSel(:,findBad) = [];

sigWidthStoreMSNSel = sigWidthStoreMSN;
sigWidthStoreMSNSel(:,findBad) = [];

%Pretty significantly different!


%Latency
%lets pull general averages, as well as BFs for FSI and MSN. This involves
%first finding the BF at each frequency for FSIs and MSNs. Use
%cleanBinFSI/MSN for this. 


latAvFSI = NaN(length(orderFile),6,2);
latFSIsingleLat = NaN(length(orderFile),6,2);
for i = 1:length(orderFile)
    %now lets find FSI BFs using cleanBinFSI
    for j = 1:threshAmp(i,2)
        [maxVal maxInd] = max(cleanBinFSI{i}(:,j));
        fsiBF(i,j+(6-threshAmp(i,2))) = maxInd;
        [maxVal maxInd] = max(cleanBinMSN{i}(:,j));
        msnBF(i,j+(6-threshAmp(i,2))) = maxInd;
    end
    combVals = cleanBinFSI{i} .* cleanBinMSN{i};
    combVals(~isnan(combVals)) = 1;
    combValStore{i} = combVals;
    areaStore(i,1) = length(find(~isnan(cleanBinFSI{i})));
    areaStore(i,2) = length(find(~isnan(cleanBinMSN{i})));
    areaStore(i,3) = length(find(~isnan(combVals)));
    
    
    %first pull FSI latmap
    tmpLatMap = latMapFSI{i};
    %make zeroes into NaNs
    tmpLatMap(tmpLatMap == 0) = NaN;
    %save average for each amplitude level
    for j = 1:threshAmp(i,2)
        latAvFSI(i,j+(6-threshAmp(i,2)),1) = nanmean(tmpLatMap(:,j));
        latAvFSI(i,j+(6-threshAmp(i,2)),2) = length(find(tmpLatMap(:,j) > 0));
        latFSIsingleLat(i,j+(6-threshAmp(i,2)),1) = tmpLatMap(fsiBF(i,j+(6-threshAmp(i,2))),j);
        latFSIsingleLat(i,j+(6-threshAmp(i,2)),2) = tmpLatMap(msnBF(i,j+(6-threshAmp(i,2))),j);
    end
    latFSIcomb{i} = tmpLatMap.*combVals;
    
    %first pull FSI latmap
    tmpLatMap = latMapMSN{i};
    %make zeroes into NaNs
    tmpLatMap(tmpLatMap == 0) = NaN;
    %save average for each amplitude level
    for j = 1:threshAmp(i,2)
        latAvMSN(i,j+(6-threshAmp(i,2)),1) = nanmean(tmpLatMap(:,j));
        latAvMSN(i,j+(6-threshAmp(i,2)),2) = length(find(tmpLatMap(:,j) > 0));
        latMSNsingleLat(i,j+(6-threshAmp(i,2)),1) = tmpLatMap(fsiBF(i,j+(6-threshAmp(i,2))),j);
        latMSNsingleLat(i,j+(6-threshAmp(i,2)),2) = tmpLatMap(msnBF(i,j+(6-threshAmp(i,2))),j);
    end
    
    latMSNcomb{i} = tmpLatMap.*combVals;
    
    subLatVal{i} = latMSNcomb{i} - latFSIcomb{i};
end

for i = 1:length(orderFile)
    valMSN = reshape(latMSNcomb{i},1,[]);
    valFSI = reshape(latFSIcomb{i},1,[]);
    try
        valLatSig(i) = signrank(valMSN,valFSI);
    catch
        valLatSig(i) = NaN;
    end
      
    meanVal(i) = nanmean(nanmean(latMSNcomb{i} - latFSIcomb{i}));
end

valLatSig(findBad) = [];
meanVal(findBad) = [];






