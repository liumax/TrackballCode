function [s] = auditoryStatsAnalysis(nexFile,fileName);

disp(fileName)

%Variables I may want to adjust

%these are variables for doing stats. 
statsWindow = [-0.500,0.500];
statsBin = [0.005, 0.01, 0.02]; %this is the bin size I want to use for statistics
statsBinNum = zeros(length(statsBin),1);
statsBinVector = cell(length(statsBin),1);
for i = 1:length(statsBin)
    statsBinNum(i) = (statsWindow(2)-statsWindow(1))/statsBin(i);
    statsBinVector{i} = [statsWindow(1)+statsBin(i)/2:statsBin(i):statsWindow(2)-statsBin(i)/2];
end

%event name detection/number of events detection
%% 
eventSize = length(nexFile.events);
eventNames = cell(eventSize,1);
eventTstamps = cell(eventSize,1);
eventTstampSize = zeros(eventSize,1);

for i= 1:eventSize
    eventNames{i} = nexFile.events{i,1}.name;
    eventTstamps{i} = nexFile.events{i,1}.timestamps;
    eventTstampSize(i) = length(eventTstamps{i});
end

%determines how many repeated events there were. This only works with
%single event type!!!
realEvents = find(eventTstampSize>1);
realEventNum = length(realEvents);


%%finding all units based on waves subfile
%This relies on having noise/multiunit be unit a
totalUnitSize = length(nexFile.waves);
unitNames = cell(totalUnitSize,1); %for storing all unit names
unitWaves = cell(totalUnitSize,1); %for storing average waveform
unitWavesSize = zeros(totalUnitSize,1); %number of points for waveform
channelID = zeros(totalUnitSize,1); %for channel number
unitLetters = zeros(totalUnitSize,1); %for numbers indicating unit identity

for i = 1:totalUnitSize
    unitNames{i} = nexFile.waves{i,1}.name;
    channelID(i) = str2num(unitNames{i}(~isletter(unitNames{i})));
    unitLetters(i) = double(unitNames{i}(end)); %this converts the letter designations to numbers
    %which start at 97 for a, increase by one for each letter
    %(b = 98, c = 99, etc)
    unitWaves{i} = nexFile.waves{i,1}.waveforms;
    unitWavesSize(i) = nexFile.waves{i,1}.NPointsWave;
end

unitNoise = find(unitLetters == 97);
unitCells = find(unitLetters ~= 97);
trueUnits = length(unitCells);

%this generates cellIndex, in which column 1 is the probe channel, column
%2 is the index of the cell in question, and column 3 is the corresponding
%noise/multiunit

cellIndex = [channelID(unitCells,1),unitCells];
multiunitIndex = [channelID(unitNoise,1),unitNoise];
for i = 1:length(cellIndex)
    cellIndex(i,3) = multiunitIndex(find(multiunitIndex(:,1) == cellIndex(i,1)),2);
end

%% finds when noise is delivered
toneTimes= nexFile.events{realEvents}.timestamps;
y=diff(toneTimes);
x=find(y>0.1);
x=x+1;
z=[toneTimes(1);toneTimes(x)];
toneTimes = z;

%%Preps cell array for firing data.
cellSize = length(nexFile.neurons);
cellData = cell(cellSize,1);

for i = 1:cellSize
    cellData{i} = nexFile.neurons{i}.timestamps;
end


%this bins data into prespecified bins

indivStatHist = cell(cellSize,length(toneTimes));%info for individual for stats

%Combines all bins from before cue, and preserves bins post cue
preCueStat = cell(cellSize,length(statsBin)); %massive dump for all bins before cue
postCueStat = cell(cellSize,length(statsBin));

%All times for histograms 
preCueTimes = cell(length(statsBin),1);
postCueTimes = cell(length(statsBin),1);

%Mean values of all bins. Precue is binned for entire pre-cue period
preCueMeans = zeros(cellSize,length(statsBin));
postCueMeans = cell(cellSize,length(statsBin));

%Difference between average firing rate per bin and average firing rate
%before tone.
postCueDiff = cell(cellSize,length(statsBin));

%responses per tone presentation, used for calculations of variance

%this code calculates the firing rate per defined bin of all statsBin
%sizes, and then splits the histogram data into the pre and post periods.
%pre data is lumped into a single large vector, while the tone presentation
%order is preserved in the post cue period. 
for k = 1:length(statsBin)
    for j = 1:cellSize
        postFiller = zeros(length(toneTimes),length(statsBinVector{k}(statsBinVector{k}>0)));
        postCueStat{j,k} = postFiller;
        for i = 1:length(toneTimes)
            %This makes all points relative to event. Also eliminates
            %points outside desired range
            toneHolder = cellData{j}(cellData{j}>toneTimes(i)+statsWindow(1) & cellData{j} < toneTimes(i) + statsWindow(2));
            toneHolder = toneHolder - toneTimes(i);
            %This will generate a histogram based on each individual trace,
            %which can be used to generate standard deviations.
            [counts,centers] = hist(toneHolder,statsBinVector{k});
            %below code necessary to prevent bugs with row vectors
            countSize = size(counts);
            centerSize = size(centers);
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            indivStatHist{j,i} = [counts'*(1/statsBin(k)),centers'];
            counts = [];
            centers = [];
            toneHolder = [];
            preCueStat{j,k} = [preCueStat{j,k};indivStatHist{j,i}(indivStatHist{j,i}(:,2)<0,1)];
            postCueStat{j,k}(i,:) = indivStatHist{j,i}(indivStatHist{j,i}(:,2)>0,1)';
        end
        preCueMeans(j,k) = mean(preCueStat{j,k});
%         clear indivStatHist
    end
    preCueTimes{k} = indivStatHist{j,i}(indivStatHist{j,i}(:,2)<0,2);
    postCueTimes{k} = indivStatHist{j,i}(indivStatHist{j,i}(:,2)>0,2);
end

%this will store alpha values for the stats bins.
binStats = cell(cellSize,length(statsBin));

%This calculates the difference of the firing rate per bin from the mean (using
%firing preceding the stimulus delivery), as well as the statistical
%difference of each individual bin from the summed bins before the cue.
for k = 1:length(statsBin)
    for j = 1:cellSize
        binStats{j,k} = zeros(length(postCueTimes{k}),1);
        postCueMeans{j,k} = zeros(length(postCueTimes{k}),1);
        postCueDiff{j,k} = zeros(length(postCueTimes{k}),1);
        for i = 1:length(postCueTimes{k})
            postCueMeans{j,k}(i) = mean(postCueStat{j,k}(:,i));
            postCueDiff{j,k}(i) = postCueMeans{j,k}(i) - preCueMeans(j,k);
            binStats{j,k}(i) = ranksum(postCueStat{j,k}(:,i),preCueStat{j,k});
        end
    end
end

    %% 
s = struct;

s.AnalysisVariables.StatsWindow = statsWindow;
s.AnalysisVariables.StatsBins = statsBin;

s.Events.EventNames = eventNames;
s.Events.EventTimeStamps = eventTstamps;
s.Events.TrialRepetitions = length(toneTimes);

s.Units.Names = unitNames;
s.Units.Waves = unitWaves;
s.Units.UnitIndex = cellIndex;

%this generates cellIndex, in which column 1 is the probe channel, column
%2 is the index of the cell in question, and column 3 is the corresponding
%noise/multiunit

s.Statistics.AlphaValues = binStats;
s.Statistics.FiringRateDifference = postCueDiff;
s.Statistics.TimePoints = postCueTimes;


end