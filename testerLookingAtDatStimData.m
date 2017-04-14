
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
% fullData.(newName) = s;
% cd(masterFolder)
% end


%pull names of recordings out
fieldNames = fields(fullData);
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
    testData = fullData.(fieldNames{i});
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
    end
    placeHolder = placeHolder + length(desigNames);
end

%lets plot some basics

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

%lets plot ratio of pre vs post vs firing rate
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
print(hFig,'lograteVsRatio1over3','-dpdf','-r0')

hFig = figure
plot(log(store.BaselineRate),store.PostBigBin' ./ store.PreThree,'k.')
hold on
plot([min(log(store.BaselineRate)) max(log(store.BaselineRate))],[1 1],'b')
title('Log Firing Rate vs Ratio of 1/3')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'lograteVsRatio1over3','-dpdf','-r0')

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

%look at post big bin and post first bin
for i = 1:placeHolder-1
    postZeroZ(i) = (store.PostZeroBin(i)-store.BaselineRate(i))/store.BaselineSTD(i);
    postBigZ(i) = (store.PostBigBin(i)-store.BaselineRate(i))/store.BaselineSTD(i);
    postThreeZ(i) =(store.PostThree(i)-store.BaselineRate(i))/store.BaselineSTD(i);
    postFourZ(i) = (store.PostFour(i)-store.BaselineRate(i))/store.BaselineSTD(i);
end

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

%combine all units
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
[B,I] = sort(store.BaselineRate);
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
    testData = fullData.(fieldNames{i});
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






