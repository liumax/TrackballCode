%This is code to analyze the batch data from my plexon recordings. This
%code should take these recordings, and output the number of cells that
%responded to the tone (based on a significance value utilizing a rank sum
%test), for those cells, the latency of the response, and how reliable that
%response was (percentage of trials with a response). 

%this should pull the number of fields in batch data structure, which is a
%proxy for the number of recordings
recordNames=fieldnames(batchDataStructure);
fieldNum = size(recordNames,1);

%placeholder for total number of all recordings (including multiunit), and
%specifically single unit recordings
totalNum = 0;
unitNum = 0;

%this should extract the number of units
for i = 1:fieldNum
    x = length(batchDataStructure.(recordNames{i}).Units.Names);
    y = size(batchDataStructure.(recordNames{i}).Units.UnitIndex,1);
    totalNum = totalNum + x;
    unitNum = unitNum + y;
end

%generating structured array to contain data properly
statsProcessing = struct('Names',[],'UnitName',[],'Multi',zeros(totalNum,1),...
    'ResponsePercentage',zeros(totalNum,1),'AlphaValues',[],'TimePoints',[],...
    'AverageFire',zeros(totalNum,1),'HistData',[]);

placeHolder = 1;

for i = 1:fieldNum
    x = length(batchDataStructure.(recordNames{i}).Units.Names);
    for j = 1:x
        statsProcessing.Names{placeHolder} = recordNames{i}; %gets name of recording
        statsProcessing.UnitName{placeHolder} =...
            batchDataStructure.(recordNames{i}).Units.Names{j};%gets name of unit
        test = strfind(statsProcessing.UnitName{placeHolder},'a');
        if test>0;
            statsProcessing.Multi(placeHolder,1) = 1;
        else
            statsProcessing.Multi(placeHolder,1) = 0;
        end
        statsProcessing.ResponsePercentage(placeHolder) = ...
            batchDataStructure.(recordNames{i}).ResponsesStats.ResponsePercentage(j);%gets percentage response
        statsProcessing.AlphaValues{placeHolder} = ...
            batchDataStructure.(recordNames{i}).Statistics.AlphaValues{j,1}; %gets alpha values at 5ms bins
        statsProcessing.TimePoints{placeHolder} = ...
            batchDataStructure.(recordNames{i}).Statistics.TimePoints{1}; %gets time points for bins
        if test > 0
            statsProcessing.AverageFire(placeHolder,1) = 5000; %avoids multiunit firing rates
        else
            statsProcessing.AverageFire(placeHolder,1) = ...
            batchDataStructure.(recordNames{i}).Plotting.AverageFiringRate(...
            find(batchDataStructure.(recordNames{i}).Units.UnitIndex(:,2) == j)); %gets average firing rates
        end
        statsProcessing.HistData{placeHolder} = ...
            batchDataStructure.(recordNames{i}).Plotting.HistogramData{j};
        placeHolder = placeHolder + 1;
    end
end

%I then want to separate out my data into multiunit and single unit data.
realUnits = struct('Names',[],'UnitName',[],...
    'ResponsePercentage',zeros(unitNum,1),...
    'AlphaValues',zeros(size(statsProcessing.AlphaValues{1},1),unitNum),...
    'TimePoints',zeros(size(statsProcessing.TimePoints{1})),...
    'AverageFire',zeros(unitNum,1),'HistData',zeros(size(statsProcessing.HistData{1},1),unitNum));

multiunits = struct('Names',[],'UnitName',[],...
    'ResponsePercentage',zeros(totalNum-unitNum,1),...
    'AlphaValues',zeros(size(statsProcessing.AlphaValues{1},1),totalNum-unitNum),...
    'TimePoints',zeros(size(statsProcessing.TimePoints{1})),...
    'HistData',zeros(size(statsProcessing.HistData{1},1),totalNum-unitNum));

realHolder = 1;
multiHolder = 1;

for i = 1:totalNum
    if statsProcessing.Multi(i) == 1;
        multiunits.Names{multiHolder} = statsProcessing.Names{i};
        multiunits.UnitName{multiHolder} = statsProcessing.UnitName{i};
        multiunits.ResponsePercentage(multiHolder) = statsProcessing.ResponsePercentage(i);
        multiunits.AlphaValues(:,multiHolder) = statsProcessing.AlphaValues{i};
        multiunits.TimePoints = statsProcessing.TimePoints{i};
        multiunits.HistData(:,multiHolder) = statsProcessing.HistData{i}(:,1);
        multiHolder = multiHolder + 1;
    elseif statsProcessing.Multi(i) == 0;
        realUnits.Names{realHolder} = statsProcessing.Names{i};
        realUnits.UnitName{realHolder} = statsProcessing.UnitName{i};
        realUnits.ResponsePercentage(realHolder) = statsProcessing.ResponsePercentage(i);
        realUnits.AlphaValues(:,realHolder) = statsProcessing.AlphaValues{i};
        realUnits.TimePoints = statsProcessing.TimePoints{i};
        realUnits.AverageFire(realHolder) = statsProcessing.AverageFire(i);
        realUnits.HistData(:,realHolder) = statsProcessing.HistData{i}(:,1);
        realHolder = realHolder + 1;
    end
end

%now we want to look at the alpha values for each unit/multiunit. The goal
%here will be to be above (or rather, below) a threshold. 
alphaHolder = realUnits.AlphaValues;
for i =1:100
    for j = 1:170
        if alphaHolder(i,j) < 0.001
            alphaHolder(i,j) = 1;
        else
            alphaHolder(i,j) = 0;
        end
    end
end

%this code finds both the first time we cross the boundary for the
%significance threshold, and the total number of time bins during which
%values are significant

diffEvents = diff(alphaHolder); %this generates diff, which picks up end of string of 1s

firstEvents = cell(170,1);
firstEnd = cell(170,1);
eventSum = zeros(170,1);


for i = 1:170
    firstEvents{i} = find(alphaHolder(:,i)==1,1,'first');
    firstEnd{i} = find(diffEvents(:,i)==-1,1,'first');
    eventSum(i) = sum(alphaHolder(:,i));
end

%this code will determine the length of all patches of 1s.
respondersLength = find(eventSum>1); %determines all the units that have more than 1 significant event

responderDur = zeros(length(respondersLength),2);

for i = 1:length(respondersLength)
    responderDur(i,1) = respondersLength(i);
    responderDur(i,2) = firstEnd{respondersLength(i)} - firstEvents{respondersLength(i)} + 1;
end

x= find(responderDur(:,2) >1);
y = find(responderDur(:,2) == 1);

for i = 1:length(y)
    figure
    plot(realUnits.HistData(:,responderDur(y(i),1)));
end

responseLags = zeros(length(x),1);
responseRel = zeros(length(x),1);
for i = 1:length(x)
    responseLags(i) = realUnits.TimePoints(firstEvents{responderDur(x(i),1)});
    responseRel(i) = realUnits.ResponsePercentage(firstEvents{responderDur(x(i),1)});
end




