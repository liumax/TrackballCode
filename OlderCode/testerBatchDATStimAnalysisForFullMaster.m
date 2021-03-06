
% %this is meant to be master code to run analysis code!
% %Establishes the home folder as the base
% masterFolder = pwd;
% masterDir = dir;
% masterDir = {masterDir.name};
% %this code here is to remove folders reported in dir like "." or ".."
% masterIndex = strfind(masterDir,'ML');
% masterIndex = find(not(cellfun('isempty', masterIndex)));
% %masterFolders is a listing of all folders that I want to target for batch
% %analysis
% masterFolders = masterDir(masterIndex);
% numFolders = size(masterFolders,2);
% for masterCount = 1:numFolders
% %generates text that signifies path to the target folder
% targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
% cd(targetPath)
% %what lists matlab files in the folder. can extract based on file type.
% test = what;
% test = test.mat;
% test = test{1};
% load(test);
% nameFind = strfind(test,'DatStimAnalysis');
% newName =strcat('x',test(1:nameFind-1));
% fullMaster.(newName) = s;
% cd(masterFolder)
% end


%pull names of recordings out
fieldNames = fields(fullMaster);
%pull individual components of names out
fieldData = cell(length(fieldNames),5);
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
end




store.FileName = cell(10,1);
store.UnitName = cell(10,1);
store.BaselineRate = zeros(10,1);
store.BaselineSTD = zeros(10,1);
store.PostZeroBin = zeros(10,1);
store.PostZeroSTD = zeros(10,1);
store.PostBigBin = zeros(10,1);
store.TrueAUC = zeros(10,1);
store.ShuffleAUC = zeros(10,2);
store.Waves = zeros(40,4,10);
placeHolder = 1;
for i = 1:length(fieldNames)
    %load data
    testData = fullMaster.(fieldNames{i});
    desigNames = testData.DesignationName;
    %pull tuning data
    %go into individual units to pull data
    for j = 1:length(desigNames)
        %pull names
        store.FileName{placeHolder - 1 + j} = fieldNames{i};
        store.UnitName{placeHolder - 1 + j} = desigNames{j};
        %pull start/end preference scores and range
        store.StartPref(placeHolder - 1 + j) = testData.RotaryData.StartPreference;
        store.EndPref(placeHolder - 1 + j) = testData.RotaryData.EndPreference;
        store.PrefRange((placeHolder - 1 + j),1) = prctile(testData.RotaryData.PreferenceDistribution,0.5);
        store.PrefRange((placeHolder - 1 + j),2) = prctile(testData.RotaryData.PreferenceDistribution,99.5);
        %pull baseline rates and std
        store.BaselineRate(placeHolder - 1 + j) = testData.(desigNames{j}).AverageRate;
        store.BaselineSTD(placeHolder - 1 + j) = testData.(desigNames{j}).AverageSTD;
        %pull spikes from selected bins
        store.PostZeroBin(placeHolder - 1 + j) = testData.(desigNames{j}).AllHistograms(81);
        store.PostZeroSTD(placeHolder - 1 + j) = testData.(desigNames{j}).HistogramStandardDeviation(81);
        store.PostBigBin(placeHolder - 1 + j) = mean(testData.(desigNames{j}).AllHistograms(81:101));
        store.PostBigBinSTD(placeHolder - 1 + j) = std(testData.(desigNames{j}).AllHistograms(81:101))';
        store.PreThree(placeHolder - 1 + j) = mean(testData.(desigNames{j}).AllHistograms(20:80));
        store.PreThreeSTD(placeHolder - 1 + j) = std(testData.(desigNames{j}).AllHistograms(20:80));
        store.PostThree(placeHolder - 1 + j) = mean(testData.(desigNames{j}).AllHistograms(81:141));
        store.PostThreeSTD(placeHolder - 1 + j) = std(testData.(desigNames{j}).AllHistograms(81:141));
        store.PostFour(placeHolder - 1 +j) = mean(testData.(desigNames{j}).AllHistograms(101:161));
        store.PostFourSTD(placeHolder - 1 +j) = std(testData.(desigNames{j}).AllHistograms(101:161));
        store.PreOne(placeHolder - 1 + j) = mean(testData.(desigNames{j}).AllHistograms(60:80));
        store.PostOne(placeHolder - 1 + j) = mean(testData.(desigNames{j}).AllHistograms(81:101));
        store.PreThreeSingle{placeHolder - 1 + j} = mean(testData.(desigNames{j}).IndividualHistograms(20:80,:));
        store.PostThreeSingle{placeHolder - 1 + j} = mean(testData.(desigNames{j}).IndividualHistograms(81:141,:));
        %pull entire histogram of response
        store.FullHist(:,placeHolder - 1 + j) = testData.(desigNames{j}).AllHistograms;
        %pull mean and std of each point along histogram. Lets re-bin into 1
        %second bins to make things cleaner and reduce the number of zeros
        for k = 1:10
            spikeSum(k,:) = sum(testData.(desigNames{j}).IndividualHistograms(1+(k-1)*20:k*20,:));
        end
        store.HistMean(:,placeHolder - 1 + j) = mean(spikeSum,2);
        store.HistVar(:,placeHolder - 1 + j) = var(spikeSum,0,2);
        store.Fano(:,placeHolder - 1 + j) = store.HistVar(:,placeHolder - 1 +j)./store.HistMean(:,placeHolder - 1 + j);
        spikeSum = [];
        %lets try and get separated examples of the interspike interval for
        %both pre and post tone periods
        negCount = 1;
        posCount = 1;
        for k = 1:size(testData.(desigNames{j}).IndividualHistograms,2)
            %find negative values first
            tarFind = find(testData.(desigNames{j}).AllRasters(:,2) == k);
            trialVals = testData.(desigNames{j}).AllRasters(tarFind,1);
            negVals = trialVals(find(trialVals < 0));
            posVals = trialVals(find(trialVals>0));
            %now need to separate out positive and negative
            if length(negVals) > 1
                negDiff(negCount:negCount + length(negVals)-2) = diff(negVals);
                negCount = negCount + length(negVals);
            end
            if length(posVals) > 1
                posDiff(posCount: posCount + length(posVals)-2) = diff(posVals);
                posCount = posCount + length(posVals);
            end
        end
        store.ISIPre{placeHolder - 1 + j} = negDiff(negDiff~=0);
        store.ISIPost{placeHolder-1 + j} = posDiff(posDiff~=0);
        negDiff = [];
        posDiff = [];
        
        if isfield(testData.(desigNames{j}),'TrueAUC')
            store.TrueAUC(placeHolder - 1 + j) = testData.(desigNames{j}).TrueAUC;
            store.ShuffleAUC(placeHolder - 1 + j,1) = prctile(testData.(desigNames{j}).ShuffleAUC,0.5);
            store.ShuffleAUC(placeHolder - 1 + j,2) = prctile(testData.(desigNames{j}).ShuffleAUC,99.5);
        end
        %store waveform
        store.Waves(:,:,placeHolder - 1 + j) = testData.(desigNames{j}).AverageWaveForms;
        
        %now lets go through and do spike shuffle
        %first need to align to lasers
        laserPulseNum = size(testData.(desigNames{j}).IndividualHistograms,2);
        %go into laser pulse data
        laserPulseTimes = testData.LaserData.LaserStartTimes;
        laserRatio = length(laserPulseTimes)/laserPulseNum;
        alignTimes = laserPulseTimes(1:laserRatio:end);
        %now take spike train and perform shuffle
%         spikeTimes = testData.(desigNames{j}).SpikeTimes;
%         spikeDiffs = diff(spikeTimes);
%         for k = 1:1000
%             randInd = randperm(length(spikeDiffs));
%             newTrain = cumsum(spikeDiffs(randInd));
%             newTrain(2:end+1) = newTrain;
%             newTrain = newTrain + spikeTimes(1);
%             %make new raster
%             newRaster = functionBasicRaster(newTrain,alignTimes,[-4 6]);
%             for l = 1:laserPulseNum
%                 %pull data from a particular trial
%                 binStore(l,1) = length(find(newRaster(:,2) == l & newRaster(:,1) < 0));
%                 binStore(l,2) = length(find(newRaster(:,2) == l & newRaster(:,1) > 0 & newRaster(:,1) < 1));
%                 binStore(l,3) = length(find(newRaster(:,2) == l & newRaster(:,1) > 1 & newRaster(:,1) < 4));
%             end
%             %now lets calculate stuff! first get baseline rate
%             baselineRateNew(placeHolder - 1 + j,k) = mean(binStore(:,1))/4;
%             baselineSTDNew(placeHolder - 1 + j,k) = std(binStore(:,1)/4);
%             postOne(placeHolder - 1 + j,k) = mean(binStore(:,2));
%             postFour(placeHolder - 1 + j,k) = mean(binStore(:,3))/3;
%         end
        placeHolder - 1 + j
            
    end
    placeHolder = placeHolder + length(desigNames);
end

[B,I] = sort(store.BaselineRate);

sortedPeaks = zeros(size(store.Waves,3),1);
sortedPeakTroughs = zeros(size(store.Waves,3),1);
sortedPeakWidth = zeros(size(store.Waves,3),1);
sortedPeakValues = zeros(size(store.Waves,3),1);
for i = 1:size(store.Waves,3);
    %open   up the set of waves
    targetWaves = squeeze(store.Waves(:,:,i));
    [maxVal maxFind] = max(max(targetWaves));
    sortedPeakValues(i) = maxVal;
    sortedWaves(:,i) = interp1([1:1:40],targetWaves(:,maxFind),[1:0.1:40],'spline');
    %turns out this seems to work fairly well. now lets pull peaks and
    %troughs
    %the peak should be the biggest point in the whole thing. Lets pull
    %that out
    [mVal,sortedPeaks(i)] = max(sortedWaves(:,i));
    %the trough should be the opposite
    [minVal,sortedPeakTroughs(i)] = min(sortedWaves(sortedPeaks(i):end,i));
    sortedPeakTroughs(i) = sortedPeakTroughs(i)/300000;
    %find peak width
    startVal = find(sortedWaves(:,i) >= maxVal/2,1,'first');
    endVal = find(sortedWaves(sortedPeaks(i):end,i) <= maxVal/2,1,'first');
    sortedPeakWidth(i) = (endVal + sortedPeaks(i) - startVal)/300000;
end

%Separate PV from MSN by rate and waveform
pvData = intersect(find(sortedPeakTroughs<0.0005),find(store.BaselineRate>2));
msnData = find(sortedPeakTroughs>0.0005);

figure
hist(sortedPeakValues,100)

figure
plot(store.BaselineRate,sortedPeakTroughs,'k.')


hFig = figure;
hold on
plot(sortedPeakWidth,sortedPeakTroughs,'k.')
% plot(sortedPeakWidth(tester),sortedPeakTroughs(tester),'r.')
% plot(sortedPeakWidth(tester2),sortedPeakTroughs(tester2),'g.')
title('Peak Half-Width vs Peak Trough')
xlabel('Peak Half-Width, seconds')
ylabel('Peak-Trough, seconds')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'PeakHalfWidthvsPeakTroughNoColor','-dpdf','-r0')

hFig = figure;
hold on
plot(store.BaselineRate,sortedPeakTroughs,'k.')
plot(store.BaselineRate(pvData),sortedPeakTroughs(pvData),'r.')
plot(store.BaselineRate(msnData),sortedPeakTroughs(msnData),'g.')
title('Firing Rate vs Peak Trough')
xlabel('Firing Rate (Hz)')
ylabel('Peak-Trough, seconds')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'BaselineRatevsPeakTrough','-dpdf','-r0')


%% Lets try separating
%make z scored averages for msns and pvs
for i = 1:length(msnData)
    msnZ(:,i) = (store.FullHist(:,msnData(i))-store.BaselineRate(msnData(i)))/ store.BaselineSTD(msnData(i));
    msnZSmooth(:,i) = smooth(msnZ(:,i),11);
end
msnZMean = mean(msnZ');
msnZSTD = std(msnZ');

hFig = figure
plot(msnZMean,'LineWidth',2)
hold on
plot(msnZMean - msnZSTD,'LineWidth',1)
plot(msnZMean + msnZSTD,'LineWidth',1)
title('Average Z-Scored Firing Rate for Putative MSNs')
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'msnZAverage','-dpdf','-r0')

for i = 1:length(pvData)
    pvZ(:,i) = (store.FullHist(:,pvData(i))-store.BaselineRate(pvData(i)))/ store.BaselineSTD(pvData(i));
    pvZSmooth(:,i) = smooth(pvZ(:,i),11);
end
pvZMean = mean(pvZ');
pvZSTD = std(pvZ');

hFig = figure
plot(pvZMean,'LineWidth',2)
hold on
plot(pvZMean - pvZSTD,'LineWidth',1)
plot(pvZMean + pvZSTD,'LineWidth',1)
title('Average Z-Scored Firing Rate for Putative PVs')
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'pvZAverage','-dpdf','-r0')

%that didnt look like terribly much. try with imagesc?
[B,Imsn] = sort(store.BaselineRate(msnData));
[B,Ipv] = sort(store.BaselineRate(pvData));


hFig = figure
imagesc(msnZSmooth(:,Imsn)',[-1 1])
colorbar
title('Z-Scored Heatmap for Putative MSNs')
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'msnZsmoothHeat','-djpeg','-r0')

hFig = figure
imagesc(pvZSmooth(:,Ipv)',[-1,1])
colorbar
title('Z-Scored Heatmap for Putative PVs')
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'pvZsmoothHeat','-djpeg','-r0')

%plot scatter of prethree vs other time points. 
statTest(1) = signrank(store.PreThree(msnData),store.PostBigBin(msnData));
statTest(2) = signrank(store.PreThree(msnData),store.PostFour(msnData));
statTest(3) = signrank(store.PostBigBin(msnData),store.PostFour(msnData));
statTest(4) = signrank(store.PreThree(pvData),store.PostBigBin(pvData));
statTest(5) = signrank(store.PreThree(pvData),store.PostFour(pvData));
statTest(6) = signrank(store.PostBigBin(pvData),store.PostFour(pvData));
statTest(7) = signrank(store.PreThree,store.PostBigBin);
statTest(8) = signrank(store.PreThree,store.PostFour);
statTest(9) = signrank(store.PostBigBin,store.PostFour);

ratioTest(1) = median(store.PreThree(msnData) ./ store.PostBigBin(msnData)');
ratioTest(2) = median(store.PreThree(msnData) ./ store.PostFour(msnData));
ratioTest(3) = median(store.PostFour(msnData) ./ store.PostBigBin(msnData)');
ratioTest(4) = median(store.PreThree(pvData) ./ store.PostBigBin(pvData)');
ratioTest(5) = median(store.PreThree(pvData) ./ store.PostFour(pvData));
ratioTest(6) = median(store.PostFour(pvData) ./ store.PostBigBin(pvData)');
ratioTest(7) = median(store.PreThree ./ store.PostBigBin');
ratioTest(8) = median(store.PreThree ./ store.PostFour);
ratioTest(9) = median(store.PostFour ./ store.PostBigBin');

ratioTestHeader = {'pre/1 msn','pre/four msn','four/1 msn','pre/1 pv','pre/4 pv','four/1 pv','pre/1tot','pre/4tot','4/1tot'};

stdTest(1) = std(store.PreThree(msnData) ./ store.PostBigBin(msnData)');
stdTest(2) = std(store.PreThree(msnData) ./ store.PostFour(msnData));
stdTest(3) = std(store.PostFour(msnData) ./ store.PostBigBin(msnData)');
stdTest(4) = std(store.PreThree(pvData) ./ store.PostBigBin(pvData)');
stdTest(5) = std(store.PreThree(pvData) ./ store.PostFour(pvData));
stdTest(6) = std(store.PostFour(pvData) ./ store.PostBigBin(pvData)');
stdTest(7) = std(store.PreThree ./ store.PostBigBin');
stdTest(8) = std(store.PreThree ./ store.PostFour);
stdTest(9) = std(store.PostFour ./ store.PostBigBin');


basicStats(1,1) = mean(store.PreThree(msnData));
basicStats(1,2) = median(store.PreThree(msnData));
basicStats(1,3) = std(store.PreThree(msnData));
basicStats(1,4) = mean(store.PostBigBin(msnData));
basicStats(1,5) = median(store.PostBigBin(msnData));
basicStats(1,6) = std(store.PostBigBin(msnData));
basicStats(1,7) = mean(store.PostFour(msnData));
basicStats(1,8) = median(store.PostFour(msnData));
basicStats(1,9) = std(store.PostFour(msnData));

basicStats(2,1) = mean(store.PreThree(pvData));
basicStats(2,2) = median(store.PreThree(pvData));
basicStats(2,3) = std(store.PreThree(pvData));
basicStats(2,4) = mean(store.PostBigBin(pvData));
basicStats(2,5) = median(store.PostBigBin(pvData));
basicStats(2,6) = std(store.PostBigBin(pvData));
basicStats(2,7) = mean(store.PostFour(pvData));
basicStats(2,8) = median(store.PostFour(pvData));
basicStats(2,9) = std(store.PostFour(pvData));

basicStats(3,1) = mean(store.PreThree);
basicStats(3,2) = median(store.PreThree);
basicStats(3,3) = std(store.PreThree);
basicStats(3,4) = mean(store.PostBigBin);
basicStats(3,5) = median(store.PostBigBin);
basicStats(3,6) = std(store.PostBigBin);
basicStats(3,7) = mean(store.PostFour);
basicStats(3,8) = median(store.PostFour);
basicStats(3,9) = std(store.PostFour);
basicStatsHeader= {'meanPre','medianPre','stdPre','mean1','median1','std1','mean14','median14','std14'};

save('statRatioTestResults.mat','statTest','ratioTest','stdTest','basicStats','basicStatsHeader','ratioTestHeader');

hFig = figure
errorbarxy(store.PreThree(msnData),store.PostFour(msnData),store.PreThreeSTD(msnData),store.PreThreeSTD(msnData),store.PostFourSTD(msnData),store.PostFourSTD(msnData),{'ko','k','k'})
hold on
plot([0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])],[0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds), Jumping 1 Putative MSNs:',num2str(statTest(2))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])])
ylim([0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostFourmsn','-dpdf','-r0')

hFig = figure
errorbarxy(store.PreThree(msnData),store.PostBigBin(msnData),store.PreThreeSTD(msnData),store.PreThreeSTD(msnData),store.PostBigBinSTD(msnData)',store.PostBigBinSTD(msnData)',{'ko','k','k'})
hold on
plot([0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])],[0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds vs 1 second) Putative MSNs:',num2str(statTest(1))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])])
ylim([0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostOnemsn','-dpdf','-r0')
%for PVs
hFig = figure
errorbarxy(store.PreThree(pvData),store.PostFour(pvData),store.PreThreeSTD(pvData),store.PreThreeSTD(pvData),store.PostFourSTD(pvData),store.PostFourSTD(pvData),{'ko','k','k'})
hold on
plot([0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])],[0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds), Jumping 1 Putative PVs:',num2str(statTest(4))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])])
ylim([0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostFourpv','-dpdf','-r0')

hFig = figure
errorbarxy(store.PreThree(pvData),store.PostBigBin(pvData),store.PreThreeSTD(pvData),store.PreThreeSTD(pvData),store.PostBigBinSTD(pvData)',store.PostBigBinSTD(pvData)',{'ko','k','k'})
hold on
plot([0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])],[0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds vs 1 second) Putative PVs:',num2str(statTest(3))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])])
ylim([0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostOnepv','-dpdf','-r0')


hFig = figure
plot(store.PreThree(msnData),store.PostFour(msnData),'k.')
hold on
plot([0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])],[0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds), Jumping 1 Putative MSNs:',num2str(statTest(2))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])])
ylim([0 max([max(store.PreThree(msnData)),max(store.PostFour(msnData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostFourmsnDOT','-dpdf','-r0')

hFig = figure
plot(store.PreThree(msnData),store.PostBigBin(msnData),'k.')
hold on
plot([0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])],[0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds vs 1 second) Putative MSNs:',num2str(statTest(1))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])])
ylim([0 max([max(store.PreThree(msnData)),max(store.PostBigBin(msnData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostOnemsnDOT','-dpdf','-r0')
%for PVs
hFig = figure
plot(store.PreThree(pvData),store.PostFour(pvData),'k.')
hold on
plot([0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])],[0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds), Jumping 1 Putative PVs:',num2str(statTest(4))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])])
ylim([0 max([max(store.PreThree(pvData)),max(store.PostFour(pvData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostFourpvDOT','-dpdf','-r0')

hFig = figure
plot(store.PreThree(pvData),store.PostBigBin(pvData),'k.')
hold on
plot([0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])],[0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])],'b')
title(strcat('Comparison of Firing Pre and Post (3 seconds vs 1 second) Putative PVs:',num2str(statTest(3))))
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])])
ylim([0 max([max(store.PreThree(pvData)),max(store.PostBigBin(pvData))])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostOnepvDOT','-dpdf','-r0')


%try plotting out scatter of z score before and after
for i = 1:length(msnData)
    zBinnedmsn(i,1) = mean(msnZ(20:80,i));
    zBinnedmsn(i,2) = mean(msnZ(81:100,i));
    zBinnedmsn(i,3) = mean(msnZ(101:160,i));
end

%this plots z score of pre vs the first second
figure
plot(zBinnedmsn(:,1),zBinnedmsn(:,2),'k.')
hold on
minVal = min([min(zBinnedmsn(:,1)),min(zBinnedmsn(:,2))]);
maxVal = max([max(zBinnedmsn(:,1)),max(zBinnedmsn(:,2))]);
plot([minVal maxVal],[minVal maxVal],'b')
xlim([minVal maxVal])
ylim([minVal maxVal])

%this plots z score of pre vs the first second
figure
plot(zBinnedmsn(:,1),zBinnedmsn(:,3),'k.')
hold on
minVal = min([min(zBinnedmsn(:,1)),min(zBinnedmsn(:,3))]);
maxVal = max([max(zBinnedmsn(:,1)),max(zBinnedmsn(:,3))]);
plot([minVal maxVal],[minVal maxVal],'b')
xlim([minVal maxVal])
ylim([minVal maxVal])

%lets do entire train shuffle, which involves some deep diving into the
%data



figure
plot(store.BaselineRate(msnData),store.PostFour(msnData) ./ store.PreThree(msnData),'k.')
hold on
plot([0 4.5],[1 1],'b')

%lets do entire train shuffle, which involves some deep diving into the
%data
%now have train shuffled data. I need to set up a loop to find what the
%percentile score is of each actual mean value relative to the distribution
%of values
for i = 1:size(postOne,1)
    whileTrig = 0;
    expVal = store.PostBigBin(i);
    shuffleArray = postOne(i,:);
    startPoint = 50;
    shuffleStep = 1;
    prevStep = 0;
    while whileTrig == 0
        testVal = prctile(shuffleArray,startPoint);
        if prevStep == 0
            if testVal > expVal
                prevStep = 1;
                startPoint = startPoint - shuffleStep;
            elseif testVal < expVal
                prevStep = -1;
                startPoint = startPoint + shuffleStep;
            elseif testVal == expVal
                postOnePrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == 1
            if testVal > expVal
                startPoint = startPoint - shuffleStep;
                if startPoint <1
                    postOnePrctile(i) = 0;
                    whileTrig = 1;
                end
            elseif testVal <= expVal
                postOnePrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == -1
            if testVal >= expVal
                postOnePrctile(i) = startPoint;
                whileTrig = 1;
            elseif testVal < expVal
                startPoint = startPoint + shuffleStep;
                if startPoint > 100
                    postOnePrctile(i) = 100;
                    whileTrig = 1;
                end
            end
        end
    end
end


for i = 1:size(baselineRateNew,1)
    whileTrig = 0;
    expVal = store.BaselineRate(i);
    shuffleArray = baselineRateNew(i,:);
    startPoint = 50;
    shuffleStep = 1;
    prevStep = 0;
    while whileTrig == 0
        testVal = prctile(shuffleArray,startPoint);
        if prevStep == 0
            if testVal > expVal
                prevStep = 1;
                startPoint = startPoint - shuffleStep;
            elseif testVal < expVal
                prevStep = -1;
                startPoint = startPoint + shuffleStep;
            elseif testVal == expVal
                BaselinePrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == 1
            if testVal > expVal
                startPoint = startPoint - shuffleStep;
                if startPoint <1
                    BaselinePrctile(i) = 0;
                    whileTrig = 1;
                end
            elseif testVal <= expVal
                BaselinePrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == -1
            if testVal >= expVal
                BaselinePrctile(i) = startPoint;
                whileTrig = 1;
            elseif testVal < expVal
                startPoint = startPoint + shuffleStep;
                if startPoint > 100
                    BaselinePrctile(i) = 100;
                    whileTrig = 1;
                end
            end
        end
    end
end


for i = 1:size(postFour,1)
    whileTrig = 0;
    expVal = store.PostFour(i);
    shuffleArray = postFour(i,:);
    startPoint = 50;
    shuffleStep = 1;
    prevStep = 0;
    while whileTrig == 0
        testVal = prctile(shuffleArray,startPoint);
        if prevStep == 0
            if testVal > expVal
                prevStep = 1;
                startPoint = startPoint - shuffleStep;
            elseif testVal < expVal
                prevStep = -1;
                startPoint = startPoint + shuffleStep;
            elseif testVal == expVal
                postFourPrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == 1
            if testVal > expVal
                startPoint = startPoint - shuffleStep;
                if startPoint <1
                    postFourPrctile(i) = 0;
                    whileTrig = 1;
                end
            elseif testVal <= expVal
                postFourPrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == -1
            if testVal >= expVal
                postFourPrctile(i) = startPoint;
                whileTrig = 1;
            elseif testVal < expVal
                startPoint = startPoint + shuffleStep;
                if startPoint > 100
                    postFourPrctile(i) = 100;
                    whileTrig = 1;
                end
            end
        end
    end
end

for i = 1:size(baselineRateNew,1)
    whileTrig = 0;
    expVal = store.PreThree(i);
    shuffleArray = baselineRateNew(i,:);
    startPoint = 50;
    shuffleStep = 1;
    prevStep = 0;
    while whileTrig == 0
        testVal = prctile(shuffleArray,startPoint);
        if prevStep == 0
            if testVal > expVal
                prevStep = 1;
                startPoint = startPoint - shuffleStep;
            elseif testVal < expVal
                prevStep = -1;
                startPoint = startPoint + shuffleStep;
            elseif testVal == expVal
                preThreePrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == 1
            if testVal > expVal
                startPoint = startPoint - shuffleStep;
                if startPoint <1
                    preThreePrctile(i) = 0;
                    whileTrig = 1;
                end
            elseif testVal <= expVal
                preThreePrctile(i) = startPoint;
                whileTrig = 1;
            end
        elseif prevStep == -1
            if testVal >= expVal
                preThreePrctile(i) = startPoint;
                whileTrig = 1;
            elseif testVal < expVal
                startPoint = startPoint + shuffleStep;
                if startPoint > 100
                    preThreePrctile(i) = 100;
                    whileTrig = 1;
                end
            end
        end
    end
end


%% lets plot some basics

hFig = figure
hist(store.BaselineRate,100);
title('Distribution of Baseline Firing Rates')

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'BaselineHistogramTerminalStim','-dpdf','-r0')


hFig = figure
errorbarxy(store.PreThree,store.PostThree,store.PreThreeSTD,store.PreThreeSTD,store.PostThreeSTD,store.PostThreeSTD,{'ko','k','k'})
hold on
plot([0 max([max(store.PreThree),max(store.PostThree)])],[0 max([max(store.PreThree),max(store.PostThree)])],'b')
title('Comparison of Firing Pre and Post (3 seconds)')
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree),max(store.PostThree)])])
ylim([0 max([max(store.PreThree),max(store.PostThree)])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostThree','-dpdf','-r0')

hFig = figure
errorbarxy(store.PreThree,store.PostFour,store.PreThreeSTD,store.PreThreeSTD,store.PostFourSTD,store.PostFourSTD,{'ko','k','k'})
hold on
plot([0 max([max(store.PreThree),max(store.PostThree)])],[0 max([max(store.PreThree),max(store.PostThree)])],'b')
title('Comparison of Firing Pre and Post (3 seconds), Jumping 1')
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree),max(store.PostFour)])])
ylim([0 max([max(store.PreThree),max(store.PostFour)])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostFour','-dpdf','-r0')

hFig = figure
errorbarxy(store.PreThree,store.PostBigBin,store.PreThreeSTD,store.PreThreeSTD,store.PostBigBinSTD',store.PostBigBinSTD',{'ko','k','k'})
hold on
plot([0 max([max(store.PreThree),max(store.PostThree)])],[0 max([max(store.PreThree),max(store.PostThree)])],'b')
title('Comparison of Firing Pre and Post (3 seconds vs 1 second)')
pbaspect([1 1 1])
xlim([0 max([max(store.PreThree),max(store.PostBigBin)])])
ylim([0 max([max(store.PreThree),max(store.PostBigBin)])])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'preThreeVSpostOne','-dpdf','-r0')

%% now lets plot out this shit as a comparison
ratioScatter(:,1) = store.PreThree ./ store.PreThree;
ratioScatter(:,2) = store.PostBigBin' ./ store.PreThree;
ratioScatter(:,3) = store.PostFour ./ store.PreThree;

%plot out as lines
hFig = figure
plot(ratioScatter')
title('Ratio of Pre vs Pre, Laser, and Post')
xlim([0.5 3.5])
set(gca,'XTick',[1:1:3])
set(gca,'XTickLabel',{'-3:0','0:1','1:4'})
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'RatioScatterLine','-dpdf','-r0')
%plot out as points
hFig = figure
plot(ratioScatter','.')
hold on
plot(mean(ratioScatter),'o')
title('Ratio of Pre vs Pre, Laser, and Post')
xlim([0.5 3.5])
set(gca,'XTick',[1:1:3])
set(gca,'XTickLabel',{'-3:0','0:1','1:4'})
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'RatioScatterPoint','-dpdf','-r0')
%plot as bar graphs

barPrep = mean(ratioScatter);
barSTD = std(ratioScatter);
barTest(1) = ranksum(ratioScatter(:,1),ratioScatter(:,2));
barTest(2) = ranksum(ratioScatter(:,1),ratioScatter(:,3));
barTest(3) = ranksum(ratioScatter(:,2),ratioScatter(:,3));

hFig = figure
bar(barPrep)
hold on
errorbarxy([1 2 3],barPrep,[0,0,0],barSTD,{'k.','k','r'})

%separate units that are suppressed vs excited in 1sec
upLaser = find(ratioScatter(:,2) > 1);
downLaser = find(ratioScatter(:,2) < 1);
%plot out separated sets
hFig = figure
plot(ratioScatter(upLaser,:)')
title('Ratio of Pre vs Pre, Laser, and Post')
xlim([0.5 3.5])
set(gca,'XTick',[1:1:3])
set(gca,'XTickLabel',{'-3:0','0:1','1:4'})
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'RatioScatterZeroPos','-dpdf','-r0')


hFig = figure
plot(ratioScatter(downLaser,:)')
title('Ratio of Pre vs Pre, Laser, and Post')
xlim([0.5 3.5])
set(gca,'XTick',[1:1:3])
set(gca,'XTickLabel',{'-3:0','0:1','1:4'})
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'RatioScatterZeroNeg','-dpdf','-r0')





%look at post big bin and post first bin
for i = 1:placeHolder-1
    postZeroZ(i) = (store.PostZeroBin(i)-store.BaselineRate(i))/store.BaselineSTD(i);
    postBigZ(i) = (store.PostBigBin(i)-store.BaselineRate(i))/store.BaselineSTD(i);
    postThreeZ(i) =(store.PostThree(i)-store.BaselineRate(i))/store.BaselineSTD(i);
    postFourZ(i) = (store.PostFour(i)-store.BaselineRate(i))/store.BaselineSTD(i);
end
%% lets plot ratio of pre vs post vs firing rate
hFig = figure
plot(store.BaselineRate,store.PostFour ./ store.PreThree,'k.')
hold on
plot([0 max(store.BaselineRate)],[1 1],'b')
title('Firing Rate vs Ratio of 4/3')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'rateVsRatio4over3','-dpdf','-r0')

hFig = figure
plot(store.BaselineRate,store.PostBigBin' ./ store.PreThree,'k.')
hold on
plot([0 max(store.BaselineRate)],[1 1],'b')
title('Firing Rate vs Ratio of 1/3')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'rateVsRatio1over3','-dpdf','-r0')

%try log of firing rate to decompress it
hFig = figure
plot(log(store.BaselineRate),store.PostFour ./ store.PreThree,'k.')
hold on
plot([min(log(store.BaselineRate)) max(log(store.BaselineRate))],[1 1],'b')
title('Log Firing Rate vs Ratio of 4/3')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'lograteVsRatio4over3','-dpdf','-r0')

hFig = figure
plot(log(store.BaselineRate),store.PostBigBin' ./ store.PreThree,'k.')
hold on
plot([min(log(store.BaselineRate)) max(log(store.BaselineRate))],[1 1],'b')
title('Log Firing Rate vs Ratio of 1/3')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'lograteVsRatio1over3','-dpdf','-r0')

%try log firing vs z score
hFig = figure
plot(log(store.BaselineRate),postFourZ,'k.')
hold on
plot([min(log(store.BaselineRate)) max(log(store.BaselineRate))],[0 0],'b')
title('Log Firing Rate vs Z 4')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'lograteVsz4','-dpdf','-r0')

hFig = figure
plot(log(store.BaselineRate),postBigZ,'k.')
hold on
plot([min(log(store.BaselineRate)) max(log(store.BaselineRate))],[0 0],'b')
title('Log Firing Rate vs Z 1')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'lograteVsz1','-dpdf','-r0')


%try log firing vs z score
hFig = figure
plot(log(store.BaselineRate(msnData)),postFourZ(msnData),'k.')
hold on
plot([min(log(store.BaselineRate)) max(log(store.BaselineRate))],[0 0],'b')
title('Log Firing Rate vs Z 4')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'lograteVsz4msn','-dpdf','-r0')

hFig = figure
plot(log(store.BaselineRate(msnData)),postBigZ(msnData),'k.')
hold on
plot([min(log(store.BaselineRate)) max(log(store.BaselineRate))],[0 0],'b')
title('Log Firing Rate vs Z 1')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'lograteVsz1msn','-dpdf','-r0')

%now try linear regression

fourThree = (store.PostFour ./ store.PreThree)';
fourThree = [ones(length(fourThree),1) fourThree];

oneThree = (store.PostBigBin' ./ store.PreThree)';
oneThree = [ones(length(oneThree),1) oneThree];

bFourThree = fourThree\log(store.BaselineRate);
bOneThree = oneThree\log(store.BaselineRate);
%I think this is getting fucked by the fact that the high firing cells skew
%everything. 

hFig = figure
plot(store.PostBigBin' ./ store.PreThree,store.PostFour ./ store.PreThree,'k.')
hold on
plot([0 max(store.PostBigBin' ./ store.PreThree)],[0 max(store.PostBigBin' ./ store.PreThree)],'b')
title('Comparison of 1/3 (x) vs 4/3 (y)')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'1v3vs4v3','-dpdf','-r0')


oneFour = (store.PostFour ./ store.PreThree)';
oneFour = [ones(length(oneFour),1) oneFour];
bOneFour = (store.PostBigBin'./store.PreThree)'\oneFour;


for i = 1:placeHolder - 1
    if store.TrueAUC(i) > store.ShuffleAUC(i,2) | store.TrueAUC(i) < store.ShuffleAUC(i,1)
        findSigAUC(i) = 1;
    else
        findSigAUC(i) = 0;
    end
end
findSigAUC = find(findSigAUC == 1);

hFig = figure
plot(store.TrueAUC,store.BaselineRate,'k.')
hold on
plot(store.TrueAUC(findSigAUC),store.BaselineRate(findSigAUC),'ro')
title({'AUC Values vs Baseline Rate';'Red Circles are Significant by 99% CI';strcat(num2str(length(findSigAUC)/(placeHolder-1)*100),'% Units Significant')})
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'BaselineVsAUCTerminal','-dpdf','-r0')


hFig = figure
hist(postZeroZ,100)
title(strcat('Histogram of Z-Score of First 50ms bin: ',num2str(median(postZeroZ)),'median'))
findBig = max(abs(postZeroZ));
xlim([-findBig findBig])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'postZeroZTerminal','-dpdf','-r0')

hFig = figure
hist(postBigZ,100)
title(strcat('Histogram of Z-Score of First 1s bin: ',num2str(median(postBigZ)),'median'))
findBig = max(abs(postBigZ));
xlim([-findBig findBig])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'postBigZTerminal','-dpdf','-r0')

hFig = figure
hist(postThreeZ,100)
title(strcat('Histogram of Z-Score of First 3s bin: ',num2str(median(postThreeZ)),'median'))
findBig = max(abs(postThreeZ));
xlim([-findBig findBig])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'postThreeZTerminal','-dpdf','-r0')

hFig = figure
hist(postFourZ,100)
title(strcat('Histogram of Z-Score of First 1-4s bin: ',num2str(median(postFourZ)),'median'))
findBig = max(abs(postFourZ));
xlim([-findBig findBig])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'postFourZTerminal','-dpdf','-r0')

%% combine all units
hFig = figure
plot(mean(store.FullHist'),'LineWidth',2)
hold on
plot(mean(store.FullHist')+std(store.FullHist')/sqrt(placeHolder-1))
plot(mean(store.FullHist')-std(store.FullHist')/sqrt(placeHolder-1))
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Averaged Firing from All Units, Aligned to Laser')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AverageAllUnitsTerminal','-dpdf','-r0')


%separate by recording
for i = 1:length(fieldNames)
    test = strfind(store.FileName,fieldNames{i});
    fieldStore = find(not(cellfun('isempty', test)));
    fieldTicker(fieldStore) = i;
end

%how about histograms of responses?
for i = 1:placeHolder - 1
    normHist(:,i) = (store.FullHist(:,i));
    zHist(:,i) = (store.FullHist(:,i)-store.BaselineRate(i))/store.BaselineSTD(i);
    smoothZ(:,i) = smooth(zHist(:,i),9);
    histMin = min(normHist(:,i));
    histMax = max(normHist(:,i));
    normHist(:,i) = (normHist(:,i)-histMin)/(histMax-histMin); 
    smoothNorm(:,i) = smooth(store.FullHist(:,i),9);
    histMin = min(smoothNorm(:,i));
    histMax = max(smoothNorm(:,i));
    smoothNorm(:,i) = (smoothNorm(:,i)-histMin)/(histMax-histMin); 
end

%for plotting separated by field/experiment
% for i = 1:length(fieldNames)
%     testData = normHist(:,fieldTicker == i);
%     testBaseline = store.BaselineRate(fieldTicker==i);
%     [B,I] = sort(testBaseline);
%     figure
%     imagesc(testData(:,I)')
% end

%order from highest baseline 

newZHist = zHist(:,I);
newNormHist = normHist(:,I);

hFig = figure
imagesc(newNormHist')
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Normalized Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapNormHistTerminal','-djpeg','-r0')

hFig = figure
imagesc(newZHist',[-2,2])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapZHistTerminal','-djpeg','-r0')


hFig = figure
imagesc(smoothNorm(:,I)')
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Normalized Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothNormHistTerminal','-djpeg','-r0')

hFig = figure
imagesc(smoothZ(:,I)',[-1,1])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score Responses Sorted by Baseline Firing Rate')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothZHistTerminal','-djpeg','-r0')
% in zscored see some pretty long lasting inhibition in the fast spiking
% cells. Could just be an artifact of the laser stimulation?

%also create alternative sorting based on spikes per 1 second after laser
[B zSort1] = sort(postBigZ);
[B zSort2] = sort(postFourZ);

[B zSort1msn] = sort(postBigZ(msnData));
[B zSort2msn] = sort(postFourZ(msnData));

[B zSort1pv] = sort(postBigZ(pvData));
[B zSort2pv] = sort(postFourZ(pvData));

hFig = figure
imagesc(smoothZ(:,msnData(zSort1msn))',[-0.5 0.5])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score MSN Responses Sorted by Response in First Sec')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothZHistmsnOrder1','-djpeg','-r0')

hFig = figure
imagesc(smoothZ(:,msnData(zSort2msn))',[-0.5 0.5])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score MSN Responses Sorted by Response in First Four')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothZHistmsnOrder2','-djpeg','-r0')

hFig = figure
imagesc(smoothZ(:,pvData(zSort1pv))',[-0.5 0.5])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score PV Responses Sorted by Response in First Sec')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothZHistpvOrder1','-djpeg','-r0')

hFig = figure
imagesc(smoothZ(:,pvData(zSort2pv))',[-0.5 0.5])
colorbar
set(gca,'XTick',[0:20:200])
set(gca,'XTickLabel',[-4:1:6])
title('Z-Score PV Responses Sorted by Response in First Four')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'HeatMapSmoothZHistpvOrder2','-djpeg','-r0')

%okay, now lets get at the running. DO we see preference for running or
%stopping?

for i = 1:placeHolder - 1
    if store.StartPref(i) > store.PrefRange(i,2)
        findSigStart(i) = 1;
    else
        findSigStart(i) = 0;
    end
end
findSigStart = find(findSigStart == 1);


for i = 1:placeHolder - 1
    if store.EndPref(i) > store.PrefRange(i,2)
        findSigEnd(i) = 1;
    else
        findSigEnd(i) = 0;
    end
end
findSigEnd = find(findSigEnd == 1);

%looks like no changes for locomotion starts or ends. What about looking at
%average locomotion trace?
for i = 1:length(fieldNames)
    testData = fullMaster.(fieldNames{i});
    %lets make average locomotion trace!
    rasterWindow = testData.Parameters.RasterWindow;
    %pull laser times
    laserTimes = testData.LaserData.LaserStartTimes;
    laserTimeDiff = min(diff(laserTimes));
    laserTimeStarts = [1;find(diff(laserTimes)>laserTimeDiff*2)+1];
    
    for j = 1:length(laserTimeStarts)
        targetTime = laserTimes(laserTimeStarts(j));
        %find start and end
        velStart = find(testData.RotaryData.Velocity(:,1)<targetTime+rasterWindow(1),1,'last');
        velEnd = velStart+((rasterWindow(2)-rasterWindow(1))/0.01);
        velSnip = testData.RotaryData.Velocity(velStart:velEnd,2);
        velHold(:,j) = velSnip;
    end
    velRaster{i} = velHold;
    velHold = [];
    laserReps(i) = length(laserTimeStarts);
    velMean(:,i) = mean(velRaster{i}');
    velSTD(:,i) = std(velRaster{i}');
    velSTE(:,i) = velSTD(:,i)/sqrt(length(laserTimeStarts));
end

%dont see much in the average velocity traces either lol


hFig = figure
hold on
for i = 1:size(velMean,2)
    plot(smooth(velMean(:,i)-mean(velMean(:,i))+(2*(i-1)),50))
end
title('Average Velocity Traces')
set(gca,'XTick',[0:100:1000])
set(gca,'XTickLabel',[-4:1:6])
xlim([0,1000])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AverageVelocityTerminal','-dpdf','-r0')

%lets try normalizing
for i = 1:size(velMean,2)
    normVel(:,i) = velMean(:,i);
    zVel(:,i) = (normVel(:,i) - mean(normVel(:,i)))/std(normVel(:,i));
    zVel(:,i) = smooth(zVel(:,i),50);
    velMax = max(normVel(:,i));
    velMin = min(normVel(:,i));
    normVel(:,i) = (normVel(:,i)-velMin)/(velMax-velMin);
    normVel(:,i) = smooth(normVel(:,i),50);
end

hFig = figure
imagesc(flipud(zVel'),[-2 2])
colorbar
title('Z-Scored Velocity')
set(gca,'XTick',[0:100:1000])
set(gca,'XTickLabel',[-4:1:6])
xlim([0,1000])
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'ZVelocitySmoothTerminal','-dpdf','-r0')

%now look at pre/post ITIs

%decide on the window I want to look at
ITIWindow = [0:0.001:0.050];
preCount = 1;
postCount = 1;

for i = 1:placeHolder -1
    %pull the data first
    preITI = store.ISIPre{i};
    postITI = store.ISIPost{i};
    fullPre(preCount:preCount + length(preITI) - 1) = preITI;
    fullPost(postCount:postCount + length(postITI) - 1) = postITI;
    
    preCount = preCount + length(preITI);
    postCount = postCount + length(postITI);
    %find number of ITIs here
    numPre(i,1) = length(preITI);
    numPost(i,1) = length(postITI);
    %remove the ITIs bigger than the ITI window
    preITI = preITI(preITI<ITIWindow(end));
    postITI = postITI(postITI<ITIWindow(end));
    %find number here
    numPre(i,2) = length(preITI);
    numPost(i,2) = length(postITI);
    
    histPre(:,i) = hist(preITI,ITIWindow);
    histPost(:,i) = hist(postITI,ITIWindow);
    
end

numPre(:,3) = numPre(:,2)./numPre(:,1);
numPost(:,3) = numPost(:,2)./numPost(:,1);

%plot overall percentage of neurons below cutoff
figure
plot(numPre(:,3),numPost(:,3),'k.')
hold on
plot([0 max(numPre(:,3))],[0 max(numPre(:,3))],'b')
title('Proportion of Spike ISIs below 50ms, Pre(x) vs Post(y)')

%plot pre vs post overall as cumulative distributions

%fullPost is max maximum, use this for limit for histogram

preCountCumDist = cumsum(hist(fullPre,[0:0.001:6]))/length(fullPre);
postCountCumDist = cumsum(hist(fullPost,[0:0.001:6]))/length(fullPost);

figure
plot(preCountCumDist)
hold on
plot(postCountCumDist,'r')
title('Cumulative Distribution of ISIs, Pre (B) vs Post (R)')

%dont seem to see much of a difference. If anything, the red line has fewer
%short ISIs and more longer ISIs.....maybe something?

%look at fano factor?
preMeanFano = mean(store.Fano(1:4,:));
postMeanFano = (store.Fano(5,:));
postPostFano = mean(store.Fano(6:end,:));
hFig = figure
plot(preMeanFano,postMeanFano,'k.')
hold on
plot([0 max(preMeanFano)],[0 max(preMeanFano)],'b')
pbaspect([1 1 1])
title('Fano Factor of Pre-Laser and 1-sec Post-Laser')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'fanoPreVs1sec','-dpdf','-r0')

hFig = figure
plot(preMeanFano,postPostFano,'k.')
hold on
plot([0 max(preMeanFano)],[0 max(preMeanFano)],'b')
pbaspect([1 1 1])
title('Fano Factor of Pre-Laser and 1-4sec Post-Laser')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'fanoPreVs1-4sec','-dpdf','-r0')






