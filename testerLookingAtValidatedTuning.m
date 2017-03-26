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

%next, go through and hunt for duplicates or overlap.
depthLim = 200;

whileTrig = 0;
searchInd = 1;
searchLim = length(fieldNames);
while whileTrig == 0
    disp(strcat('Looking at File #',num2str(searchInd)))
    %reset search limits
    searchLim = length(fieldNames);
    
    if searchInd == searchLim - 1;
        whileTrig = 1;
    end
    %first see if dates are the same
    dateComp = strcmp(fieldData{searchInd,1},fieldData{searchInd+1,1});
    if dateComp == 1
        %then compare mouse ID
        nameComp = strcmp(fieldData{searchInd,2},fieldData{searchInd+1,2});
        if nameComp == 1
            %then compare side
            sideComp = strcmp(fieldData{searchInd,3},fieldData{searchInd+1,3});
            if sideComp == 1
                %compare depths
                
                depth1 = str2double(fieldData{searchInd,4});
                depth2 = str2double(fieldData{searchInd+1,4});
                depthDiff = abs(depth1 - depth2);
                if depthDiff < depthLim
                    %initiate text to eliminate one of the options
                    disp('Overlap Detected with Files')
                    disp(strcat(fieldNames{searchInd},' or ',fieldNames{searchInd+1}))
                    whileCounter = 0; %this is the counter that gets updated to exit the loop

                    while whileCounter == 0
                        try
                            prompt = 'Decide on Units: k: keep both f: keep first s: keep second';
                            str = input(prompt,'s');
                            if str~='k' & str~='f' & str~='s'
                                error
                            else
                                whileCounter = 1;
                            end
                        catch
                        end
                    end
                    if strfind(str,'k')
                        disp('Keep Both Units')
                        searchInd = searchInd + 1;
                    elseif strfind(str,'f')
                        disp('Keeping First, Deleting Second')
                        %delete the row below search index
                        fieldData(searchInd + 1,:) = [];
                        fieldNames(searchInd + 1) = [];
                        %advance searchInd
                        searchInd = searchInd + 1;
                    elseif strfind(str,'s')
                        disp('Keeping Second, Deleting First')
                        %delete the row below search index
                        fieldData(searchInd,:) = [];
                        fieldNames(searchInd) = [];
                        %advance searchInd
                        searchInd = searchInd + 1;
                    end
                else
                    searchInd = searchInd + 1;
                end
            else
                searchInd = searchInd + 1;
            end
        else
            searchInd = searchInd + 1;
        end
    else
        searchInd = searchInd + 1;
    end
end

%Now, we should have an annotated list of the files we want to pull from
%fullMaster. 

%using these, lets go through a couple of big for loops and extract relevant details

for i = 1:length(fieldNames)
    disp(strcat('Analyzing Locomotion For ',fieldNames{i}))
    %load data from the correct field
    testData = fullMaster.(fieldNames{i});
    desigNames = testData.DesignationName;
    %first, check if there is locomotor data for ROC analysis, and see if
    %this analysis has been performed. If not, execute the ROC analysis
    locoCheck = length(testData.RotaryData.LocoStarts);
    if locoCheck > 0
        disp('Locomotion Data Found, checking for ROC analysis')
        rocCheck = isfield(testData.(desigNames{1}),'TrueAUC');
        if rocCheck == 0
            disp('Performing ROC analysis of Locomotion/Firing Rate')
            for j = 1:length(desigNames)
                [testData] = functionLocomotionROC(testData,desigNames{j});
                fullMaster.(fieldNames{i}) = testData;
            end
            
        elseif rocCheck == 1
            disp('ROC Data Detected, no further analysis needed')
        end
    elseif locoCheck == 0
        disp('No locomotor data, not checking for ROC')
    end
end
disp('Locomotor Repairs Finished')

%next, lets try and determine whether there is finer tuning related
%information in the files. if it doesnt exist, then execute finer tuning.

for i = 1:length(fieldNames)
    disp(strcat('Analyzing Tuning Responses For ',fieldNames{i}))
    %load data
    testData = fullMaster.(fieldNames{i});
    desigNames = testData.DesignationName;
    %check if there is analysis for type of tuning
    tuneCheck = isfield(testData.(desigNames{1}),'TuningType');
    if tuneCheck == 0;
        disp('Tuning Type Analysis Not Performed, Initiating')
        [decisionTuning,tuningType] = functionTuningSelectionTool(testData,fieldNames{i});
        fullMaster.(fieldNames{i}).TuningDecision = decisionTuning;
        fullMaster.(fieldNames{i}).TuningType = tuningType;
    else
        disp('Tuning Type Analysis Already Done')
    end
end
disp('Tuning Repairs Finished')

%next, lets go through the dataset and pull out information we care about

store.FileName = cell(100,1);
store.UnitName = cell(100,1);
store.DecisionTuning = zeros(100,1);
store.TuningType = cell(100,1);
store.BaselineRate = zeros(100,1);
store.TrueAUC = zeros(100,1);
store.ShuffleAUC = zeros(100,2);
store.MinLat = zeros(100,1);

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
        
    end

    placeHolder = placeHolder + length(desigNames);
end


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

