%This code is meant to look at validated tuning analysis



cd E:\170219ValidatedTuningAnalysis\Processed\reSorted
load('FullTuningAnalysis07-Mar-2017.mat')

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
% 
% %next, go through and hunt for duplicates or overlap.
% depthLim = 200;
% 
% whileTrig = 0;
% searchInd = 1;
% searchLim = length(fieldNames);
% while whileTrig == 0
%     disp(strcat('Looking at File #',num2str(searchInd)))
%     %reset search limits
%     searchLim = length(fieldNames);
%     
%     if searchInd == searchLim - 1;
%         whileTrig = 1;
%     end
%     %first see if dates are the same
%     dateComp = strcmp(fieldData{searchInd,1},fieldData{searchInd+1,1});
%     if dateComp == 1
%         %then compare mouse ID
%         nameComp = strcmp(fieldData{searchInd,2},fieldData{searchInd+1,2});
%         if nameComp == 1
%             %then compare side
%             sideComp = strcmp(fieldData{searchInd,3},fieldData{searchInd+1,3});
%             if sideComp == 1
%                 %compare depths
%                 
%                 depth1 = str2double(fieldData{searchInd,4});
%                 depth2 = str2double(fieldData{searchInd+1,4});
%                 depthDiff = abs(depth1 - depth2);
%                 if depthDiff < depthLim
%                     %initiate text to eliminate one of the options
%                     disp('Overlap Detected with Files')
%                     disp(strcat(fieldNames{searchInd},' or ',fieldNames{searchInd+1}))
%                     whileCounter = 0; %this is the counter that gets updated to exit the loop
% 
%                     while whileCounter == 0
%                         try
%                             prompt = 'Decide on Units: k: keep both f: keep first s: keep second';
%                             str = input(prompt,'s');
%                             if str~='k' & str~='f' & str~='s'
%                                 error
%                             else
%                                 whileCounter = 1;
%                             end
%                         catch
%                         end
%                     end
%                     if strfind(str,'k')
%                         disp('Keep Both Units')
%                         searchInd = searchInd + 1;
%                     elseif strfind(str,'f')
%                         disp('Keeping First, Deleting Second')
%                         %delete the row below search index
%                         fieldData(searchInd + 1,:) = [];
%                         fieldNames(searchInd + 1) = [];
%                         %advance searchInd
%                         searchInd = searchInd + 1;
%                     elseif strfind(str,'s')
%                         disp('Keeping Second, Deleting First')
%                         %delete the row below search index
%                         fieldData(searchInd,:) = [];
%                         fieldNames(searchInd) = [];
%                         %advance searchInd
%                         searchInd = searchInd + 1;
%                     end
%                 else
%                     searchInd = searchInd + 1;
%                 end
%             else
%                 searchInd = searchInd + 1;
%             end
%         else
%             searchInd = searchInd + 1;
%         end
%     else
%         searchInd = searchInd + 1;
%     end
% end
% 
% %Now, we should have an annotated list of the files we want to pull from
% %fullMaster. 
% 
% %using these, lets go through a couple of big for loops and extract relevant details
% 
% for i = 1:length(fieldNames)
%     disp(strcat('Analyzing Locomotion For ',fieldNames{i}))
%     %load data from the correct field
%     testData = fullMaster.(fieldNames{i});
%     desigNames = testData.DesignationName;
%     %first, check if there is locomotor data for ROC analysis, and see if
%     %this analysis has been performed. If not, execute the ROC analysis
%     locoCheck = length(testData.RotaryData.LocoStarts);
%     if locoCheck > 0
%         disp('Locomotion Data Found, checking for ROC analysis')
%         rocCheck = isfield(testData.(desigNames{1}),'TrueAUC');
%         if rocCheck == 0
%             disp('Performing ROC analysis of Locomotion/Firing Rate')
%             for j = 1:length(desigNames)
%                 [testData] = functionLocomotionROC(testData,desigNames{j});
%                 fullMaster.(fieldNames{i}) = testData;
%             end
%             
%         elseif rocCheck == 1
%             disp('ROC Data Detected, no further analysis needed')
%         end
%     elseif locoCheck == 0
%         disp('No locomotor data, not checking for ROC')
%     end
% end
% disp('Locomotor Repairs Finished')
% 
% %next, lets try and determine whether there is finer tuning related
% %information in the files. if it doesnt exist, then execute finer tuning.
% 
% for i = 1:length(fieldNames)
%     disp(strcat('Analyzing Tuning Responses For ',fieldNames{i}))
%     %load data
%     testData = fullMaster.(fieldNames{i});
%     desigNames = testData.DesignationName;
%     %check if there is analysis for type of tuning
%     tuneCheck = isfield(testData,'TuningType');
%     if tuneCheck == 0;
%         disp('Tuning Type Analysis Not Performed, Initiating')
%         [decisionTuning,tuningType] = functionTuningSelectionTool(testData,fieldNames{i});
%         fullMaster.(fieldNames{i}).TuningDecision = decisionTuning;
%         fullMaster.(fieldNames{i}).TuningType = tuningType;
%     else
%         disp('Tuning Type Analysis Already Done')
%     end
% end
% disp('Tuning Repairs Finished')

%next, lets go through the dataset and pull out information we care about

store.FileName = cell(100,1);
store.UnitName = cell(100,1);
store.DecisionTuning = zeros(100,1);
store.TuningType = cell(100,1);
store.BaselineRate = zeros(100,1);
store.TrueAUC = zeros(100,1);
store.ShuffleAUC = zeros(100,2);
store.MinLat = zeros(100,1);
store.Waves = zeros(40,4,100);
% store.LatMaps = cell;

placeHolder = 1;

for i = 1:length(fieldNames)
    %load data
    testData = fullMaster.(fieldNames{i});
    desigNames = testData.DesignationName;
    %pull tuning data
    store.DecisionTuning(placeHolder:placeHolder + length(desigNames)-1) = testData.TuningDecision;
    store.TuningType(placeHolder:placeHolder + length(desigNames)-1) = testData.TuningType;
    %go into individual units to pull data
    for j = 1:length(desigNames)
        store.FileName{placeHolder - 1 + j} = fieldNames{i};
        store.UnitName{placeHolder - 1 + j} = desigNames{j};
        store.BaselineRate(placeHolder - 1 + j) = testData.(desigNames{j}).AverageRate;
        if length(min(testData.(desigNames{j}).LatencyMap(testData.(desigNames{j}).LatencyMap>0))) > 0
            store.MinLat(placeHolder - 1 + j) = min(testData.(desigNames{j}).LatencyMap(testData.(desigNames{j}).LatencyMap>0));
        end
        
        if isfield(testData.(desigNames{j}),'TrueAUC')
            store.TrueAUC(placeHolder - 1 + j) = testData.(desigNames{j}).TrueAUC;
            store.ShuffleAUC(placeHolder - 1 + j,1) = prctile(testData.(desigNames{j}).ShuffleAUC,0.5);
            store.ShuffleAUC(placeHolder - 1 + j,2) = prctile(testData.(desigNames{j}).ShuffleAUC,99.5);
        end
        store.LatMaps{placeHolder - 1 + j} = testData.(desigNames{j}).LatencyMap;
        %store waveform
        store.Waves(:,:,placeHolder - 1 + j) = testData.(desigNames{j}).AverageWaveForms;
        
    end

    placeHolder = placeHolder + length(desigNames);
end

%lets go through the waveforms and try and pick out representative ones

%as a first pass, lets go through all the shits and select the peak wave as
%the model
% sortedWaves = zeros(40,size(store.Waves,3));
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


tester = intersect(find(sortedPeakTroughs<0.00046),find(sortedPeakWidth<0.00014));
tester2 = intersect(find(sortedPeakTroughs<0.00046),find(sortedPeakWidth<0.0002 & sortedPeakWidth>0.00014));
tester3 = [1:1:length(sortedPeakWidth)];
selUnique = unique([tester;tester2]);
tester3(selUnique) = [];

figure
hist(sortedPeakValues,100)

figure
plot(store.BaselineRate,sortedPeakTroughs,'k.')


hFig = figure;
hold on
plot(sortedPeakWidth,sortedPeakTroughs,'k.')
plot(sortedPeakWidth(tester),sortedPeakTroughs(tester),'r.')
plot(sortedPeakWidth(tester2),sortedPeakTroughs(tester2),'g.')
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
plot(store.BaselineRate(tester),sortedPeakTroughs(tester),'r.')
plot(store.BaselineRate(tester2),sortedPeakTroughs(tester2),'g.')
title('Firing Rate vs Peak Trough')
xlabel('Firing Rate (Hz)')
ylabel('Peak-Trough, seconds')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'BaselineRatevsPeakTrough','-dpdf','-r0')

hFig = figure
hold on
stem3(sortedPeakWidth(tester3),sortedPeakTroughs(tester3),store.BaselineRate(tester3),'k','filled')
stem3(sortedPeakWidth(tester),sortedPeakTroughs(tester),store.BaselineRate(tester),'r','filled')
stem3(sortedPeakWidth(tester2),sortedPeakTroughs(tester2),store.BaselineRate(tester2),'g','filled')
xlabel('FWHM (sec)')
ylabel('Peak-Trough (sec)')
zlabel('Baseline Rate (Hz)')
grid on

%% lets try and pull out latencies

latStore = zeros(100,1);
counter = 1;
for i = 1:length(store.FileName)
    testMap = store.LatMaps{i};
    findLat = find(testMap);
    latVals = testMap(findLat);
    latStore(counter:counter+length(latVals)-1) = latVals;
    counter = counter + length(latVals);
end

figure
hist(latStore(latStore<0.15),[0:0.001:0.15])

%lets try only for tuned neurons
latStoreSel = zeros(100,1);
counter = 1;
tuneFind = find(store.DecisionTuning);
for i = 1:length(tuneFind)
    testMap = store.LatMaps{tuneFind(i)};
    findLat = find(testMap);
    latVals = testMap(findLat);
    latStoreSel(counter:counter+length(latVals)-1) = latVals;
    counter = counter + length(latVals);
end

figure
hist(latStoreSel(latStoreSel<0.15),[0:0.001:0.15])

%This overall looks fairly reasonable. I'm seeing most things with above 6
%ms latency, and some shit with short latency, which is somewhat expected
%based on how I calculate things. 

%Now lets go through the data and get the tuning information from each unit
%so I can get a sense of tuning bandwidths.

placeHolder = 1;
latStore = struct;
for i = 1:length(fieldNames)
    %load data
    testData = fullMaster.(fieldNames{i});
    desigNames = testData.DesignationName;
    %pull tuning data
    DecisionTuning = find(testData.TuningDecision);
    %go into individual tuned units to pull data
    for j = 1:length(DecisionTuning)
        latStore(placeHolder).LatMap = testData.(desigNames{DecisionTuning(j)}).LatencyMap;
        latStore(placeHolder).Freqs = testData.SoundData.UniqueFrequencies;
        latStore(placeHolder).DBs = testData.SoundData.UniqueDBs;
        latStore(placeHolder).MaxAmp = find(testData.(desigNames{DecisionTuning(j)}).LatencyMap(:,end));
        placeHolder = placeHolder + 1;
    end
end

%% Testing extraction of responses to given frequencies


%now pull all units and frequencies to get sense of how many cells were
%presented with specific frequencies
placeHolder = 1;
for i = 1:length(fieldNames)
    %load data
    testData = fullMaster.(fieldNames{i});
    desigNames = testData.DesignationName;
    %pull tuning data
    
    %go into individual tuned units to pull data
    for j = 1:length(desigNames)
        freqStore(placeHolder:placeHolder - 1 + length(testData.SoundData.UniqueFrequencies)) = testData.SoundData.UniqueFrequencies;
        placeHolder = placeHolder + length(testData.SoundData.UniqueFrequencies);
    end
end

[fullFreqSet ia ic] = unique(freqStore);
%find abundance of each thing!
for i = 1:length(fullFreqSet)
    fullFreqAb(i) = length(find(ic == i));
end

%turns out there is a bit of a bug with the frequencies. Round, eliminate,
%and group
fullFreqSet = round(fullFreqSet);
[fullFreqSet ia ic] = unique(fullFreqSet);
newFreqAb = zeros(length(fullFreqSet),1);
for i = 1:length(ic)
    newFreqAb(ic(i)) = newFreqAb(ic(i))+ fullFreqAb(i);
end


%lets just get readout of frequencies represented with responses at peak
%amp
counter = 1;
for i = 1:length(latStore)
    targetFreqs = latStore(i).Freqs(latStore(i).MaxAmp);
    if length(targetFreqs > 0)
        respFreqs(counter:counter + length(targetFreqs) - 1) = targetFreqs;
        counter = counter + length(targetFreqs);
    end
end

[respFreqSet ia ic] = unique(respFreqs);

for i = 1:length(respFreqSet)
    respFreqAb(i) = length(find(ic == i));
end

%turns out there is a bit of a bug with the frequencies. Round, eliminate,
%and group
respFreqSet = round(respFreqSet);
[respFreqSet ia ic] = unique(respFreqSet);
newRespFreqAb = zeros(length(respFreqSet),1);
for i = 1:length(ic)
    newRespFreqAb(ic(i)) = newRespFreqAb(ic(i))+ respFreqAb(i);
end

respRatio = newRespFreqAb ./ newFreqAb; %this plots the ratio of the responses vs presentations relative to number of cells

figure
semilogx(respFreqSet,respRatio)

%% this more or less works! Lets try to do it in a for loop to figure things
%out. 

%okay, what we are going to do is as follows: pull all frequencies, pull
%all decibels, and then determine the total range. This will allow me to
%generate a big array of m x n, with decibels and frequencies. Then I can
%go through the dataset and pull individual sets and use those to fill one
%array for occurrence and one array for responses
freqCounter = 1;
dbCounter = 1;
% allFreqs = zeros(100,1);
% allDBS = zeros(100,1);
for i = 1:length(latStore)
    targetFreqs = latStore(i).Freqs;
    targetDBs = latStore(i).DBs;
    allFreqs(freqCounter:freqCounter + length(targetFreqs) - 1) = targetFreqs;
    allDBs(dbCounter:dbCounter + length(targetDBs) - 1) = targetDBs;
    freqCounter = freqCounter + length(targetFreqs);
    dbCounter = dbCounter + length(targetDBs);
end

allFreqs = unique(round(allFreqs));
allDBs = unique(round(allDBs));

%based on these, we can generate an array of dbs x freqs

incidenceArray = zeros(length(allDBs),length(allFreqs));
responseArray = zeros(length(allDBs),length(allFreqs));

for i = 1:length(latStore)
    %pull freqs and dbs
    targetFreqs = round(latStore(i).Freqs);
    targetDBs = round(latStore(i).DBs);
    %match to overall lists
    [Cf,ia,ib] = intersect(allFreqs,targetFreqs);
    [Cd,ic,id] = intersect(allDBs,targetDBs);
    %fill incidence array.
    for j = 1:length(ic)
        for k = 1:length(ia)
            incidenceArray(ic(j),ia(k)) = incidenceArray(ic(j),ia(k)) + 1;
        end
    end
    incidenceArray(ic,ia) = incidenceArray(ic,ia) + 1;
    %now need to fill out response array
    latFind = find(latStore(i).LatMap);
    %now I need to find the coordinates
    coordFreq = rem(latFind,length(targetFreqs));
    coordFreq(coordFreq == 0) = length(targetFreqs);
    coordDB = fix((latFind-1)/length(targetFreqs))+1;
    coordDB(coordDB<1) = 1;
    for j = 1:length(coordDB)
        responseArray(ic(coordDB(j)),ia(coordFreq(j))) = responseArray(ic(coordDB(j)),ia(coordFreq(j))) + 1;
    end
end

perRespArray = responseArray ./ incidenceArray;

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
imagesc(perRespArray)
colormap('parula')
colorbar
set(gca,'XTick',[1:1:17])
set(gca,'XTickLabel',allFreqs)
set(gca,'YTick',[1:1:6])
set(gca,'YTickLabel',allDBs)
title('Heatmap Of Percent Responses Relative To Amp/Freq')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'heatmapResponsePercentage','-dpdf','-r0')

%% Now that we've done that, lets try something a bit more complex: lets use the same code to pull out the relevant response latencies

incidenceArray = zeros(length(allDBs),length(allFreqs));
responseArray = ones(length(allDBs),length(allFreqs));
responseLatStore = cell(length(allDBs),length(allFreqs));
for i = 1:length(latStore)
    %pull freqs and dbs
    targetFreqs = round(latStore(i).Freqs);
    targetDBs = round(latStore(i).DBs);
    %match to overall lists
    [Cf,ia,ib] = intersect(allFreqs,targetFreqs);
    [Cd,ic,id] = intersect(allDBs,targetDBs);
    %fill incidence array.
    for j = 1:length(ic)
        for k = 1:length(ia)
            incidenceArray(ic(j),ia(k)) = incidenceArray(ic(j),ia(k)) + 1;
        end
    end
    incidenceArray(ic,ia) = incidenceArray(ic,ia) + 1;
    %now need to fill out response array
    latFind = find(latStore(i).LatMap);
    %now I need to find the coordinates
    coordFreq = rem(latFind,length(targetFreqs));
    coordFreq(coordFreq == 0) = length(targetFreqs);
    coordDB = fix((latFind-1)/length(targetFreqs))+1;
    coordDB(coordDB<1) = 1;
    for j = 1:length(coordDB)
        responseLatStore{ic(coordDB(j)),ia(coordFreq(j))}(responseArray(ic(coordDB(j)),ia(coordFreq(j)))) = latStore(i).LatMap(latFind(j));
        responseArray(ic(coordDB(j)),ia(coordFreq(j))) = responseArray(ic(coordDB(j)),ia(coordFreq(j))) + 1;
    end
end

%lets pull the averages
responseLatAverage = zeros(length(allDBs),length(allFreqs));
for i = 1:length(allDBs)
    for j = 1:length(allFreqs)
        responseLatAverage(i,j) = mean(responseLatStore{i,j});
    end
end

figure
imagesc(responseLatAverage,[0 0.05])
colormap('parula')
colorbar

%% plot out average waveforms

testWaves = store.Waves(:,:,tester);
testWaves2 = store.Waves(:,:,tester2);
testWaves3 = store.Waves(:,:,tester3);
% 
% figure
% plot(squeeze(mean(testWaves,2)))
% figure
% plot(squeeze(mean(testWaves2,2)))
% figure
% plot(squeeze(mean(testWaves3,2)))

figure
hold on
plot(squeeze(mean(squeeze(mean(testWaves,2)),2)),'r','LineWidth',2)
plot(squeeze(mean(squeeze(mean(testWaves2,2)),2)),'g','LineWidth',2)
plot(squeeze(mean(squeeze(mean(testWaves3,2)),2)),'k','LineWidth',2)
title('Average Waveform Per Grouping (PV (R),UNK (G), MSN (K)')

%lets look at mean rates
baseRateUntuned = store.BaselineRate(store.DecisionTuning == 0);
baseRateTuned = store.BaselineRate(store.DecisionTuning == 1);

histBaseUn = hist(baseRateUntuned,[0:0.5:25]);
histBaseTu = hist(baseRateTuned,[0:0.5:25]);

figure
plot([0:0.5:25],histBaseUn,[0:0.5:25],histBaseTu)
%look pretty similar

meanUntuned = mean(baseRateUntuned);
meanTuned = mean(baseRateTuned);

%how about locomotion tuning?
store.AUCSig = zeros(length(store.FileName),1);
for i = 1:length(store.FileName);
    if store.TrueAUC(i) > store.ShuffleAUC(i,2) | store.TrueAUC(i) < store.ShuffleAUC(i,1)
        store.AUCSig(i) = 1;
    end
    
end

%find where there is an AUC score calculated at all
store.AUCFind = zeros(280,1);
store.AUCFind(store.TrueAUC>0) = 1;

%pull this together
aucTuningMatrix = zeros(sum(store.AUCFind),4);
aucTuningMatrix(:,1) = store.DecisionTuning(store.AUCFind == 1);
aucTuningMatrix(:,2) = store.TrueAUC(store.AUCFind == 1);
aucTuningMatrix(:,3) = store.AUCSig(store.AUCFind == 1);
aucTuningMatrix(:,4) = store.BaselineRate(store.AUCFind == 1);

%see relationship of AUC scores vs tuning
aucSigUn = aucTuningMatrix(aucTuningMatrix(:,1) == 0,3);
aucSigTu = aucTuningMatrix(aucTuningMatrix(:,1) == 1,3);

%looks like more significant cells in tuned cells? Also may be biased by
%crappy untuned units

figure
plot(aucTuningMatrix(:,2),aucTuningMatrix(:,4),'b.');
hold on
plot(aucTuningMatrix(aucTuningMatrix(:,3) == 1,2),aucTuningMatrix(aucTuningMatrix(:,3) == 1,4),'r.');

%dont see anything super consistent there...


%what about using this as a dataset to pull up good parameters for tuning?

%lets look at number of consecutive bin and peak rate. Best would be to
%divide by frequency. 
tunedData = struct;
findTuned = find(store.DecisionTuning == 1);

for i = 1:length(findTuned)
    %basically want to go through the full dataset and pull the relevant
    %stuff from it
    
    %basically, it seems like I was doing analysis based on 5%, 1%, and
    %0.1% significance. Will follow the 1% line for excitation, 5% line for
    %inhibition
   
    %first look at overall histogram
    tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).OverallHist = fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).AllHistogramSig.Histogram;
    %pull mean value
    tunedData.MeanRate(i) = fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).AllHistogramSig.MeanBaseline;
    %pull peak value
    tunedData.OverallMax(i) = max(max(tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).OverallHist(:,1)));
    %pull min value during histogram period
    tunedData.OverallMin(i) = min(min(tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).OverallHist(:,1)));
    %pull longest consecutive significant peak duration
    testHold = tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).OverallHist(:,3);
    testDiff = diff(testHold);
    testDiffPos = find(testDiff == 1);
    testDiffNeg = find(testDiff == -1);
    %need to see if first value is significant. if so, then need to insert
    %a 0 value for the first value of testDiffPos
    if testHold(1) == 1
        testDiffPos(2:end+1) = testDiffPos;
        testDiffPos(1) = 0;
    end
    
    if ~isempty(testDiffPos) & ~isempty(testDiffNeg)
        if length(testDiffPos) == length(testDiffNeg)
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        elseif length(testDiffPos) > length(testDiffNeg)
            
            testDiffPos = testDiffPos(1:length(testDiffNeg));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        else length(testDiffPos) < length(testDiffNeg)

            testDiffNeg = testDiffNeg(1:length(testDiffPos));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        end
    else
        testDur = 0;
    end
    tunedData.PosDuration(i) = testDur;
    %now do for negative. Go with 5% line
    testHold = tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).OverallHist(:,5);
    testDiff = diff(testHold);
    testDiffPos = find(testDiff == 1);
    testDiffNeg = find(testDiff == -1);
    %need to see if first value is significant. if so, then need to insert
    %a 0 value for the first value of testDiffPos
    if testHold(1) == 1
        testDiffPos(2:end+1) = testDiffPos;
        testDiffPos(1) = 0;
    end
    
    if ~isempty(testDiffPos) & ~isempty(testDiffNeg)
        if length(testDiffPos) == length(testDiffNeg)
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        elseif length(testDiffPos) > length(testDiffNeg)
            testDiffPos = testDiffPos(1:length(testDiffNeg));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        else length(testDiffPos) < length(testDiffNeg)
            testDiffNeg = testDiffNeg(1:length(testDiffPos));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        end
    else
        testDur = 0;
    end
    tunedData.NegDuration(i) = testDur;
    
    tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).SigBins = fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).SpecHistogramSig;
    %need to get into each of these values. 
    tester = reshape(fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).SpecHistogramSig,[],1);
    for j = 1:length(tester)
        
    end
    
    tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).PeakResp = fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).PeakMap;
    tunedData.MaxPeak(i) = max(max(tunedData.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).PeakResp));
    %find histogram zero and tone end
    histZero = find(fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).HistBinVector > 0,1,'first');
    histToneEnd = find(fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).HistBinVector < fullMaster.(store.FileName{findTuned(i)}).SoundData.ToneDuration,1,'last');
    tunedData.MinVal(i) = min(min(min(fullMaster.(store.FileName{findTuned(i)}).(store.UnitName{findTuned(i)}).FreqDBHistograms(:,:,histZero:histToneEnd))));
end

%find all excited cells
exciteFind = strfind(store.TuningType,'e');
Index = find(not(cellfun('isempty', exciteFind)));
%match this to the indices among tuned neurons
[sharedVals,exciteInd] = intersect(find(store.DecisionTuning == 1),Index);
%here idxsIntoA is the correct values for excited neurons

%find all excited cells
inhibFind = strfind(store.TuningType,'i');
Index = find(not(cellfun('isempty', inhibFind)));
%match this to the indices among tuned neurons
[sharedVals,inhibInd] = intersect(find(store.DecisionTuning == 1),Index);
%here idxsIntoA is the correct values for excited neurons


figure
plot(tunedData.OverallMax./tunedData.MeanRate,tunedData.PosDuration,'k.')
hold on
plot(tunedData.OverallMax(exciteInd)./tunedData.MeanRate(exciteInd),tunedData.PosDuration(exciteInd),'r.')
plot(tunedData.OverallMax(inhibInd)./tunedData.MeanRate(inhibInd),tunedData.PosDuration(inhibInd),'g.')

figure
plot(tunedData.OverallMax-tunedData.MeanRate,tunedData.PosDuration,'k.')
hold on
plot(tunedData.OverallMax(exciteInd)-tunedData.MeanRate(exciteInd),tunedData.PosDuration(exciteInd),'r.')
plot(tunedData.OverallMax(inhibInd)-tunedData.MeanRate(inhibInd),tunedData.PosDuration(inhibInd),'g.')

%now find non-responsive cells
%lets look at number of consecutive bin and peak rate. Best would be to
%divide by frequency. 
untunedData = struct;
findUntuned = find(store.DecisionTuning == 0);

for i = 1:length(findUntuned)
    %basically want to go through the full dataset and pull the relevant
    %stuff from it
    
    %basically, it seems like I was doing analysis based on 5%, 1%, and
    %0.1% significance. Will follow the 1% line for excitation, 5% line for
    %inhibition
   
    %first look at overall histogram
    untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).OverallHist = fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).AllHistogramSig.Histogram;
    %pull mean value
    untunedData.MeanRate(i) = fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).AllHistogramSig.MeanBaseline;
    %pull peak value
    untunedData.OverallMax(i) = max(max(untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).OverallHist(:,1)));
    %pull min value during histogram period
    untunedData.OverallMin(i) = min(min(untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).OverallHist(:,1)));
    %pull longest consecutive significant peak duration
    testHold = untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).OverallHist(:,3);
    testDiff = diff(testHold);
    testDiffPos = find(testDiff == 1);
    testDiffNeg = find(testDiff == -1);
    %need to see if first value is significant. if so, then need to insert
    %a 0 value for the first value of testDiffPos
    if testHold(1) == 1
        testDiffPos(2:end+1) = testDiffPos;
        testDiffPos(1) = 0;
    end
    
    if ~isempty(testDiffPos) & ~isempty(testDiffNeg)
        if length(testDiffPos) == length(testDiffNeg)
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        elseif length(testDiffPos) > length(testDiffNeg)
            testDiffPos = testDiffPos(1:length(testDiffNeg));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        else length(testDiffPos) < length(testDiffNeg)
            testDiffNeg = testDiffNeg(1:length(testDiffPos));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        end
    else
        testDur = 0;
    end
    untunedData.PosDuration(i) = testDur;
    %now do for negative. Go with 5% line
    testHold = untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).OverallHist(:,5);
    testDiff = diff(testHold);
    testDiffPos = find(testDiff == 1);
    testDiffNeg = find(testDiff == -1);
    %need to see if first value is significant. if so, then need to insert
    %a 0 value for the first value of testDiffPos
    if testHold(1) == 1
        testDiffPos(2:end+1) = testDiffPos;
        testDiffPos(1) = 0;
    end
    
    if ~isempty(testDiffPos) & ~isempty(testDiffNeg)
        if length(testDiffPos) == length(testDiffNeg)
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        elseif length(testDiffPos) > length(testDiffNeg)
            testDiffPos = testDiffPos(1:length(testDiffNeg));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        else length(testDiffPos) < length(testDiffNeg)
            testDiffNeg = testDiffNeg(1:length(testDiffPos));
            testDiffPos = reshape(testDiffPos,[],1);
            testDiffNeg = reshape(testDiffNeg,[],1);
            testDur = testDiffNeg - testDiffPos;
            testDur = max(testDur);
        end
    else
        testDur = 0;
    end
    untunedData.NegDuration(i) = testDur;
    
    untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).SigBins = fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).SpecHistogramSig;
    %need to get into each of these values. 
    tester = reshape(fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).SpecHistogramSig,[],1);
    for j = 1:length(tester)
        
    end
    
    untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).PeakResp = fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).PeakMap;
    untunedData.MaxPeak(i) = max(max(untunedData.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).PeakResp));
    %find histogram zero and tone end
    histZero = find(fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).HistBinVector > 0,1,'first');
    histToneEnd = find(fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).HistBinVector < fullMaster.(store.FileName{findUntuned(i)}).SoundData.ToneDuration,1,'last');
    untunedData.MinVal(i) = min(min(min(fullMaster.(store.FileName{findUntuned(i)}).(store.UnitName{findUntuned(i)}).FreqDBHistograms(:,:,histZero:histToneEnd))));
end




