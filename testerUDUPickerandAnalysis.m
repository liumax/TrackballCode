%This is meant to be code to analyze the s struct file left from the basic
%analysis. It is meant to isolate only the tuned cells and perform analyses
%on them

%placeholder for opening file
fileName = '';

%pull tuned cells
indexTuned = find(s.decisionTuning == 1);

%pull names of tuned cells
namesTuned = s.DesignationName(indexTuned);

%generate structured array for storage of data
newS = struct;
newS.DesignationName = namesTuned;

%move to new structured array
for i = 1:length(namesTuned)
    newS.(namesTuned{i}) = s.(namesTuned{i});
end

%compute whether spikes match up or not.

numUnits = length(namesTuned); %pulls number of units
numCross = nchoosek(numUnits,2); %pulls # of combinations 
crossDesig = nchoosek([1:1:numUnits],2); %figures out all possible matches
zeroVals = zeros(numCross,3); %generates matrix for storing data

for i = 1:numCross
    x = s.(s.DesignationName{crossDesig(i,1)}).SpikeTimes; %pulls spike times from one unit
    y = s.(s.DesignationName{crossDesig(i,2)}).SpikeTimes; %pulls spike times from another unit
    
    if ~isempty(length(intersect(x,y)))
        zeroVals(i,1) = length(intersect(x,y)); %saves number of shared spikes
        zeroVals(i,2) = zeroVals(i,1)/length(x); %computes percentage of x spikes that are shared
        zeroVals(i,3) = zeroVals(i,1)/length(y); %computes percentage of y spikes that are shared
    end
    
end

%need to build out something to analyze/match up responses?






