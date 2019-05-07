%This code is janky code to look at identified pairs. The goal is the
%produce plots of the cross correlograms (without smoothing),
%representations of responses, and things like difference in best
%frequency.

%% Extraction

targetFiles = {'180315_ML180306C_R17_3218mid1_fullTuningAnalysis.mat',...
    '180718_ML180619B_L_AudStr_pen1_3000_fullTuningAnalysis.mat',...
    '180718_ML180619C_R_AudStr_pen1_3000_fullTuningAnalysis.mat',...
    '190123ML181105E_RAudStr3526pen1rec1tuningAndDMRAnalysis.mat',...
    '190205ML181105C_RAudStr3667pen2rec1tuningDMRAnalysis.mat',...
    '190206ML181105F_RAudStr3633pen2rec1tuningDMRAnalysis.mat',...%After this is new stuff added in!
    '180614_ML180515A_R17_2886_fullTuningAnalysis.mat',...
    '180622_ML180515C_R17_pen3_2640_fullTuningAnalysis.mat',...
    '180622_ML180515D_R17_3300_secondfullTuningAnalysis.mat',...
    '180717_ML180619A_R_AudStr_pen2_2850_fullTuningAnalysis.mat',...
    '190405ML190307B_RAudStrPen1Rec1_3526TuningDMRAnalysis.mat',...
    '190409ML190307E_RAudStrPen1Rec1_3570TuningDMRAnalysis.mat',...
    '190418ML190307G_RAudStrPen1Rec2_3546TuningDMRAnalysis.mat',...
    '190419ML190307J_RAudStrPen1Rec1_3585TuningDMRAnalysis.mat'};


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
6	283	95 %after this is new data
7   7  283 %180614
8 9 3 %180622 180515C
9 218 8 %180622 180515D
9 218 22
9 218 68
9 218 74
10 16 265%180717
11 274 276
12 271 239
12 271 267
12 271 274
13 280 273
13 280 287
13 280 298
14 0 82];

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
    numFreqs = s.SoundData.NumFreqs;
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
%         %grab trial by trial data so I can do that kind of noise
%         %correlation analysis
%         tmpDataFSI = [];
%         tmpDataMSN = [];
%         for k = 1:numDBs
%             for l = 1+s.SoundData.WhiteNoise:numFreqs
%                 %load the correct data
%                 tmpDataFSI(:,k,l-s.SoundData.WhiteNoise) = s.(tester{fsiInd(j)}).LatPeakBin{l,k}.BinnedSpikesTone-s.(tester{fsiInd(j)}).LatPeakBin{l,k}.BinnedSpikesToneBase;
%                 tmpDataMSN(:,k,l-s.SoundData.WhiteNoise) = s.(tester{msnInd(j)}).LatPeakBin{l,k}.BinnedSpikesTone-s.(tester{msnInd(j)}).LatPeakBin{l,k}.BinnedSpikesToneBase;
%             end
%         end
%         trialDataFSI{overCount + j - 1} = tmpDataFSI;
%         trialDataMSN{overCount + j - 1} = tmpDataMSN;
        
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
    shank1order = sort(numName(s.PeakChanVals <= 32));
%     if i == 1
%         shank1order = sort(numName(s.PeakChanVals <= 32));
%     else
%         shank1order = sort(numName(s.PeakChanVals < 32));
%     end
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

%% Now lets process. Start by trying to find good/bad

%now lets find things!
counter = 1;
findBad = [];
minVal = 5;
for i = 1
    tmpVals = sigNum(:,i);
    tmpVals = find(tmpVals < minVal);
    findBad(counter:counter+length(tmpVals)-1) = tmpVals;
    counter = counter + length(tmpVals);
end

findBad = unique(findBad);

%make good index
findGood = [1:length(orderFile)];
findGood(findBad) = [];

%manual annotation
manBad = [4,6,9,11,12,13,38];
manGood = [1:length(orderFile)];
manGood(manBad) = [];

%% Calculate asymmetry based on published calculations

%now lets calculate asymmetry indices! This is A = (R-L)/(R+L), with R
%being area under right hand side, and L being area under L side.

Rs = sum(truexCorr(11:20,:));
Ls = sum(truexCorr(22:31,:));

asymmetryVal = (Rs - Ls)./(Rs + Ls);

Rs = (truexCorr(11:20,:));
Ls = (truexCorr(22:31,:));

asymmetryVal = (Rs - Ls)./(Rs + Ls);
asymmetryVal = abs(asymmetryVal);
asymmetryVal = mean(asymmetryVal);
    
%it looks like the asymmetry index is at best a poor method to look at this
%kind of issue. 

%% Try and plot things out. We want to figure out a thing to separate what I think are good vs bad units.

%just plot everything
%plot out true correlations
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
    subplot(borderControl,borderControl,i)
    hold on
    plot(xAxVals,truexCorr(:,i),'k','LineWidth',2)
    
    plot(xAxVals,rangeBounds(:,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,2,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,1,i),'Color',[0.7 0.7 0.7])
    title(num2str(asymmetryVal(i)))
end


spikeGraphName = 'AllCrossCorrSig';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%plot out true correlations with labeling with color. 
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
    subplot(borderControl,borderControl,i)
    hold on
    if ismember(i,manGood)
        plot(xAxVals,truexCorr(:,i),'r','LineWidth',2)
    else
        plot(xAxVals,truexCorr(:,i),'k','LineWidth',2)
    end
    plot(xAxVals,rangeBounds(:,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,2,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,1,i),'Color',[0.7 0.7 0.7])
    title(num2str(asymmetryVal(i)))
end


spikeGraphName = 'AllCrossCorrSigLabel';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

% %lets reflect over zero to see symmetry. Keep the post-zero in the correct
% %direction, reverse the other side. 
% %plot out true correlations
% hFig = figure;
% set(hFig, 'Position', [10 80 1800 1200])
% % borderControl = ceil(sqrt(length(orderFile)));
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
% xAxVals = [-0.01:0.0005:0.01];
% for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
% %     subplot(6,4,i)
%     hold on
%     plot(xAxVals(21:end),truexCorr(21:end,i),'r','LineWidth',2)
%     plot(xAxVals(21:end),truexCorr(21:-1:1,i),'k','LineWidth',2)
%     if ismember(i,manGood)
%         title(num2str(asymmetryVal(i)),'Color','r')
%     else
%         title(num2str(asymmetryVal(i)))
%     end
%     
% end
% 
% 
% spikeGraphName = 'AllCrossCorrSigReflect';
% savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')


%okay, lets try subtracting the shuffled baseline and see what things look
%like.

subCorr = truexCorr - rangeBounds;
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
    subplot(borderControl,borderControl,i)
%     subplot(6,4,i)
    hold on
    plot(xAxVals,subCorr(:,i),'k','LineWidth',2)
    if ismember(i,manGood)
        title(num2str(asymmetryVal(i)),'Color','r')
    else
        title(num2str(asymmetryVal(i)))
    end
    
end


spikeGraphName = 'AllCrossCorrBaseSub';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


% %Try smoothing by 3
% smoothWindCorr = 3;
% hFig = figure;
% set(hFig, 'Position', [10 80 1800 1200])
% % borderControl = ceil(sqrt(length(orderFile)));
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
% xAxVals = [-0.01:0.0005:0.01];
% for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
% %     subplot(6,4,i)
%     hold on
%     plot(xAxVals,smooth(subCorr(:,i),smoothWindCorr),'k','LineWidth',2)
%     if ismember(i,manGood)
%         title(num2str(asymmetryVal(i)),'Color','r')
%     else
%         title(num2str(asymmetryVal(i)))
%     end
%     
% end
% 
% 
% spikeGraphName = 'AllCrossCorrBaseSubSmooth3';
% savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')
% 
% %try doing reflection of smoothed by 3, baseline subtracted
% hFig = figure;
% set(hFig, 'Position', [10 80 1800 1200])
% % borderControl = ceil(sqrt(length(orderFile)));
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
% xAxVals = [-0.01:0.0005:0.01];
% for i = 1:length(orderFile);
%     subplot(borderControl,borderControl,i)
% %     subplot(6,4,i)
%     tester = smooth(subCorr(:,i),smoothWindCorr);
%     hold on
%     plot(xAxVals(21:end),tester(21:end),'r','LineWidth',2)
%     plot(xAxVals(21:end),tester(21:-1:1),'k','LineWidth',2)
%     if ismember(i,manGood)
%         title(num2str(asymmetryVal(i)),'Color','r')
%     else
%         title(num2str(asymmetryVal(i)))
%     end
%     
% end
% 
% 
% spikeGraphName = 'AllCrossCorrBaseSubSmooth3SigReflect';
% savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')

%Okay, maybe solution is this (minus baseline, smooth three, then look for
%point by point differences). I need to align to whatever peak there is
%though. 

hasPeak = [3,4,6,9,11,12,13,15,17,19,20,22,23,24,28,30,33,38,39];
smoothWindCorr = 3;
for i = 1:length(orderFile)
    smoothBaseSubxCorr(:,i) = smooth(subCorr(:,i),smoothWindCorr);
end

[pks findPeaks] = max(smoothBaseSubxCorr);

%try doing reflection around peak (if there) of smoothed by 3, baseline subtracted
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
    subplot(borderControl,borderControl,i)
%     subplot(6,4,i)
    tester = smooth(subCorr(:,i),smoothWindCorr);
    hold on
    if ismember(i,hasPeak);
        plot(xAxVals(findPeaks(i):findPeaks(i) + 10),tester(findPeaks(i):findPeaks(i)+10),'r.-','LineWidth',2)
        plot(xAxVals(findPeaks(i):findPeaks(i) + 10),tester(findPeaks(i):-1:findPeaks(i)-10),'k','LineWidth',2)
        plot(xAxVals(findPeaks(i):findPeaks(i) + 10),tester(findPeaks(i):-1:findPeaks(i)-10),'g*','LineWidth',2)
    else
        plot(xAxVals(21:31),tester(21:31),'r','LineWidth',2)
        plot(xAxVals(21:31),tester(21:-1:11),'k','LineWidth',2)
    end
    
    if ismember(i,manGood)
        title(num2str(asymmetryVal(i)),'Color','r')
    else
        title(num2str(asymmetryVal(i)))
    end
    
end

%I think some of these arent doing proper peak alignment. Manually adjust. 

peakAdj = zeros(length(orderFile),1);
peakAdj(9) = -0.5;
peakAdj(11) = -0.5;
peakAdj(12) = +1;
peakAdj(13) = +0.5;
peakAdj(22) = +0.5;
peakAdj(23) = +1;
peakAdj(30) = -0.5;

alertPeak = find(abs(peakAdj) == 0.5);

dist = 10;
for i = 1:length(orderFile)
    tester = smoothBaseSubxCorr(:,i);
    if ismember(i,hasPeak) && ismember(i,alertPeak)
        alignPoint = ceil(findPeaks(i) + peakAdj(i));
        reflectCorr(:,i,1) = tester(alignPoint:alignPoint+dist);
        reflectCorr(:,i,2) = tester(alignPoint-1:-1:alignPoint-1-dist);
    elseif ismember(i,hasPeak) && ~ismember(i,alertPeak)
        alignPoint = findPeaks(i) + peakAdj(i);
        reflectCorr(:,i,1) = tester(alignPoint:alignPoint+dist);
        reflectCorr(:,i,2) = tester(alignPoint:-1:alignPoint-dist);
    else
        reflectCorr(:,i,1) = tester(21:21+dist);
        reflectCorr(:,i,2) = tester(21:-1:21-dist);
    end
end


%extract data from after 1ms and up to 4 ms out. This I can then average
%etc. 
reflectCorrSelVal = reflectCorr(3:7,:,:);
reflectCorrSelValMean = squeeze(mean(reflectCorrSelVal));
%problem: negative values cause problems for calculation. 
reflectCorrSelValMean(:,1) = reflectCorrSelValMean(:,1) + mean(rangeBounds)';
reflectCorrSelValMean(:,2) = reflectCorrSelValMean(:,2) + mean(rangeBounds)';
newAsymm = (reflectCorrSelValMean(:,1) - reflectCorrSelValMean(:,2))./(reflectCorrSelValMean(:,1) + reflectCorrSelValMean(:,2));


hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
    subplot(borderControl,borderControl,i)
%     subplot(6,4,i)
    tester = smooth(subCorr(:,i),smoothWindCorr);
    hold on
    if ismember(i,hasPeak);
        plot(xAxVals(21:21+dist),reflectCorr(:,i,1),'r.-','LineWidth',2)
        plot(xAxVals(21:21+dist),reflectCorr(:,i,2),'k','LineWidth',2)
%         plot(xAxVals(21:21+dist),reflectCorr(:,i,2),'g*','LineWidth',2)
    else
        plot(xAxVals(21:21+dist),reflectCorr(:,i,1),'r','LineWidth',2)
        plot(xAxVals(21:21+dist),reflectCorr(:,i,2),'k','LineWidth',2)
    end
    
    if ismember(i,manGood)
        title(num2str(newAsymm(i)),'Color','r')
    else
        title(num2str(newAsymm(i)))
    end
    
end
spikeGraphName = 'AllCrossProperAlignedFoldedNoGreen';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
    subplot(borderControl,borderControl,i)
%     subplot(6,4,i)
    tester = smooth(subCorr(:,i),smoothWindCorr);
    hold on
    if ismember(i,hasPeak);
        plot(xAxVals(21:21+dist),reflectCorr(:,i,1),'r.-','LineWidth',2)
        plot(xAxVals(21:21+dist),reflectCorr(:,i,2),'k','LineWidth',2)
        plot(xAxVals(21:21+dist),reflectCorr(:,i,2),'g*','LineWidth',2)
    else
        plot(xAxVals(21:21+dist),reflectCorr(:,i,1),'r','LineWidth',2)
        plot(xAxVals(21:21+dist),reflectCorr(:,i,2),'k','LineWidth',2)
    end
    
    if ismember(i,manGood)
        title(num2str(newAsymm(i)),'Color','r')
    else
        title(num2str(newAsymm(i)))
    end
    
end
spikeGraphName = 'AllCrossProperAlignedFolded';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

% This is not completely good. I think the current asymmetry index is
% getting closer but not quite there. 



%plot out tone vs non-tone
hFig = figure;
set(hFig, 'Position', [10 80 1800 1200])
% borderControl = ceil(sqrt(length(orderFile)));
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
xAxVals = [-0.01:0.0005:0.01];
for i = 1:length(orderFile);
    subplot(borderControl,borderControl,i)
%     subplot(6,4,i)
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

% 
% %plot out correlations with synchronous spikes removed
% hFig = figure;
% set(hFig, 'Position', [10 80 1800 1200])
% % borderControl = ceil(sqrt(length(orderFile)));
% subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
% xAxVals = [-0.01:0.0005:0.01];
% for i = 1:length(orderFile);
% %     subplot(borderControl,borderControl,i)
%     subplot(6,4,i)
%     hold on
%     plot(xAxVals,remStore(:,i),'k','LineWidth',2)
%     plot(xAxVals,rangeBounds(:,i),'Color',[0.7 0.7 0.7])
%     plot(xAxVals,rangeBoundsPRCTL(:,2,i),'Color',[0.7 0.7 0.7])
%     plot(xAxVals,rangeBoundsPRCTL(:,1,i),'Color',[0.7 0.7 0.7])
% end
% 
% 
% spikeGraphName = 'AllCrossCorrSigSynchRemoved';
% savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')



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

tarLimsX = [-0.1 0.2];
for i = 1:length(orderFile)
    hFig = figure;
    set(hFig, 'Position', [10 80 1500 1000])
%     borderControl = ceil(sqrt(length(orderFile)));
    subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
    subplot(3,3,1)
    hold on
    if ismember(i,manGood)
        plot(xAxVals,truexCorr(:,i),'r','LineWidth',2)
    else
        plot(xAxVals,truexCorr(:,i),'k','LineWidth',2)
    end
    plot(xAxVals,rangeBounds(:,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,2,i),'Color',[0.7 0.7 0.7])
    plot(xAxVals,rangeBoundsPRCTL(:,1,i),'Color',[0.7 0.7 0.7])
    
    %plot overlaid rasters
    subplot(2,3,4)
    hold on
    plot(rastersFSI{i}(:,1),rastersFSI{i}(:,3),'r.')
    plot(rastersMSN{i}(:,1),rastersMSN{i}(:,3),'k.')
    xlim(tarLimsX)
    ylim([0 max(rastersFSI{i}(:,3))])
    
    %now plot rasters
    subplot(2,3,2)
    plot(rastersFSI{i}(:,1),rastersFSI{i}(:,3),'k.')
    xlim(tarLimsX)
    ylim([0 max(rastersFSI{i}(:,3))])
    title('Raster FSI')
    subplot(2,3,5)
    plot(rastersMSN{i}(:,1),rastersMSN{i}(:,3),'k.')
    title('Raster MSN')
    xlim(tarLimsX)
    ylim([0 max(rastersFSI{i}(:,3))])
    subplot(1,3,3)
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
    
    
    spikeGraphName = strcat('RastFrequencyHist-',num2str(i));
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
end
% 
% 
% %Plot trial by trial correlations?
% sigCutoff = 0.05;
% 
% for i = 1:length(orderFile)
%     counter = 1;
%     counter2 = 1;
%     numFreqs = size(trialDataFSI{i},3);
%     numDBs = size(trialDataMSN{i},2);
%     fsiData = trialDataFSI{i};
%     msnData = trialDataMSN{i};
%     hFig = figure;
%     set(hFig, 'Position', [10 10 1500 800])
% %     borderControl = ceil(sqrt(length(orderFile)));
%     subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.05], [0.05 0.05]);
%     for j = 1:numFreqs
%         for k = 1:numDBs
%             subplot(numDBs,numFreqs,(k-1)*numFreqs + j)
%             hold on
%             %store all spikes
%             bigStore(counter:counter + length(fsiData(:,k,j)) - 1,:) = [fsiData(:,k,j)-mean(fsiData(:,k,j)),msnData(:,k,j)-mean(msnData(:,k,j))];
%             counter = counter + length(fsiData(:,k,j));
%             %store just significant FSI response spikes
%             sigTar = signStore{i,1};
%             
%             if sigTar(j,k) < sigCutoff
%                 bigStore2(counter2:counter2 + length(fsiData(:,k,j)) - 1,:) = [fsiData(:,k,j)-mean(fsiData(:,k,j)),msnData(:,k,j)-mean(msnData(:,k,j))];
%                 counter2 = counter2 + length(fsiData(:,k,j));
%             end
%             
%             
%             plot(fsiData(:,k,j),msnData(:,k,j),'k.')
%             [b,bintr,bintjm] = gmregress(fsiData(:,k,j),msnData(:,k,j),sigCutoff);
%             plot([min(fsiData(:,k,j)) max(fsiData(:,k,j))],[min(fsiData(:,k,j))*b(2) + b(1) max(fsiData(:,k,j))*b(2) + b(1)],'r')
%             set(gca,'xtick',[])
%             set(gca,'xticklabel',[])
%             set(gca,'ytick',[])
%             set(gca,'yticklabel',[])
%         end
%     end
%     
%     spikeGraphName = strcat('TrialCorr-',num2str(i));
%     savefig(hFig,spikeGraphName);
% 
%     %save as PDF with correct name
%     set(hFig,'Units','Inches');
%     pos = get(hFig,'Position');
%     set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(hFig,spikeGraphName,'-dpdf','-r0')
%     
%     %keep bigstore data.
%     storeBigStore{i} = bigStore;
%     storeBigStore2{i} = bigStore2;
%     
% end
% 
% for i = 1:24
%     tester = corrcoef(storeBigStore{i}(:,1),storeBigStore{i}(:,2));
%     corrTrial(i) = tester(2);
% end
% for i = 1:24
%     tester = corrcoef(storeBigStore2{i}(:,1),storeBigStore2{i}(:,2));
%     corrTrial2(i) = tester(2);
% end



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
sigWidthStoreSel(:,manBad) = [];

sigWidthStoreMSNSel = sigWidthStoreMSN;
sigWidthStoreMSNSel(:,manBad) = [];

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
    %lets just look at 70dB
    valMSN = latMSNcomb{i}(:,end);
    valFSI = latFSIcomb{i}(:,end);
    try
        valLatSig70dB(i) = signrank(valMSN,valFSI);
    catch
        valLatSig70dB(i) = NaN;
    end
      
    meanVal70dB(i) = nanmean(nanmean(latMSNcomb{i} - latFSIcomb{i}));
end

valLatSig(findBad) = [];
meanVal(findBad) = [];


%% Looking at trial by trial differences.
%for this analysis, we will look at latPeakBin data. 






