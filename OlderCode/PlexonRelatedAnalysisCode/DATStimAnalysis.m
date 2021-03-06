function [s] = DATStimAnalysis(fileName);
% 
% fileName = 'ML150730B_AP14_DV2913_whitenoise'

%Reads in NEX File
% [nexFile] = readNexFile(strcat(char(fileName),'.nex'));
[nexFile] = readNexFile(fileName);

%Variables I may want to adjust
rasterWindow = [-2,4];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
histBin = 0.25; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)-histBin/2:histBin:rasterWindow(2)-histBin/2];

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

%find start time
x = strfind(eventNames,'Start');
y = find(not(cellfun('isempty', x)));
startTime = eventTstamps{y};

x = [];
y = [];

%find end time
x = strfind(eventNames,'Stop');
y = find(not(cellfun('isempty', x)));
stopTime = eventTstamps{y};

x = [];
y = [];

totalTime = stopTime - startTime; %total recording time in seconds. 
%precision down to 1/10th of millisecond

%determines how many repeated events there were
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
for i = 1:size(cellIndex,1)
    cellIndex(i,3) = multiunitIndex(find(multiunitIndex(:,1) == cellIndex(i,1)),2);
end

%%
%This is to compute average firing rate of every 'true' unit.
averageFiringRate = zeros(trueUnits,1);
for i = 1:trueUnits
    averageFiringRate(i) = length(nexFile.neurons{cellIndex(i,2),1}.timestamps)/totalTime;
end


%% 

%% 
toneTimes= nexFile.events{realEvents}.timestamps;
y=diff(toneTimes);
x=find(y>0.1);
x=x+1;
z=[toneTimes(1);toneTimes(x)];
toneTimes = z;

%Preps cell array for firing data.
cellSize = length(nexFile.neurons);
cellData = cell(cellSize,1);

for i = 1:cellSize
    cellData{i} = nexFile.neurons{i}.timestamps;
end

%Generates raster plot for tone times as well as histogram
rasterHolder = 1;
toneRaster = zeros(100000,2);
masterToneRaster = cell(cellSize,1);
masterToneHist = cell(cellSize,1);
indivToneHist = cell(cellSize,length(toneTimes));
%cycles through cells
for j = 1:cellSize
    %cycles through tones presentation times. extracts info for
    %individual and overall rasters.
    for i = 1:length(toneTimes)
        %This makes all points relative to event. Also eliminates
        %points outside desired range
        toneHolder = cellData{j}(cellData{j}>toneTimes(i)+rasterWindow(1) & cellData{j} < toneTimes(i) + rasterWindow(2));
        toneHolder = toneHolder - toneTimes(i);
        %This will generate a histogram based on each individual trace,
        %which can be used to generate standard deviations.
        [counts,centers] = hist(toneHolder,histBinVector);
        %below code necessary to prevent bugs with row vectors
        countSize = size(counts);
        centerSize = size(centers);
        if countSize(1)>countSize(2)
            counts = counts';
        end
        if centerSize(1)>centerSize(2)
            centers = centers';
        end
        indivToneHist{j,i} = [counts'*(1/histBin),centers'];
        counts = [];
        centers = [];
        %fills in large raster plot. holder updated position.
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
        rasterHolder = rasterHolder + length(toneHolder);
        toneRaster(toneRaster(:,1) == 0,:) = [];
        toneHolder = [];
    end
    masterToneRaster{j} = toneRaster;
    toneRaster = zeros(100000,2);
    [counts,centers] = hist(masterToneRaster{j}(:,2),histBinVector);
    countSize = size(counts);
    centerSize = size(centers);
    if countSize(1)>countSize(2)
        counts = counts';
    end
    if centerSize(1)>centerSize(2)
        centers = centers';
    end
    masterToneHist{j} = [counts'*(1/histBin)/length(toneTimes),centers'];
    counts = [];
    centers = [];
    
       
end

 %here, I need to basically z-score the data, then make sure all changes are
%aligned in the same direction

%this finds average firing rate for the entire period before stim (in this
%case, 2 seconds), and uses that to calculate average rates and standard
%deviation. This then zscores the firing rate by subtracting mean, dividing
%by std.
for j = 1:size(cellIndex,1)
    zScoreFire = zeros(length(histBinVector),size(cellIndex,1));
    for i = 1:size(cellIndex,1)
        x = find(histBinVector<0);%find period before laser stim
        aveFire = mean(masterToneHist{cellIndex(i,2)}(x,1));
        aveStd = std(masterToneHist{cellIndex(i,2)}(x,1));
        zScoreFire = (masterToneHist{cellIndex(i,2)}(:,1)-aveFire)/aveStd;
        zScoreAve = mean(zScoreFire(find(histBinVector >0 & histBinVector<2)));
    end
    
    holder = struct;
    holder.AnalysisWindow = rasterWindow;
    holder.AnalysisBins = histBin;

    holder.Events.EventNames = eventNames;
    holder.Events.EventTimeStamps = eventTstamps;
    holder.Events.Trials = length(toneTimes);

    holder.Units.Names = unitNames;
    holder.Units.AverageWaves = unitWaves;
    holder.Units.AverageFiringRate = averageFiringRate;
    holder.Units.UnitIndex = cellIndex;

    holder.Plotting.Rasters = masterToneRaster;
    holder.Plotting.Histograms = masterToneHist;
    holder.Plotting.ZScoreFiring = zScoreFire;
    holder.Plotting.ZScoreAverage = zScoreAve;
    
    s{j} = holder;
end

stdHolder = zeros(length(histBinVector),length(toneTimes));
steHolder = zeros(length(histBinVector),cellSize);

for i = 1:cellSize
    for j = 1:length(toneTimes)
        stdHolder(:,j) = indivToneHist{i,j}(:,1);
    end
    steHolder(:,i) = std(stdHolder,0,2)/sqrt(length(toneTimes));
end

stePlotter = zeros(length(histBinVector),cellSize,2);
for i = 1:cellSize
    stePlotter(:,i,1) = masterToneHist{i,1}(:,1)-steHolder(:,i);
    stePlotter(:,i,2) = masterToneHist{i,1}(:,1)+steHolder(:,i);
end

set(0, 'DefaulttextInterpreter', 'none')

for i = 1:size(cellIndex,1)
    hFig = figure;
    set(hFig,'Units','inches');
    set(hFig,'Position',[1 1 6 8]);

    subplot(3,2,1)
    plot(unitWaves{cellIndex(i,2)})
    set(gca, 'Units', 'inches');
    set(gca,'OuterPosition',[0 5.5 3 2.7]);
    title('Average Waveform')

    mTextBox = uicontrol('style','text');
    descr = {'File:';
        fileName;
        'Unit:';
        unitNames{cellIndex(i,2)}
        'Average Firing Rate (Hz):';
        averageFiringRate(i)};
    set(mTextBox,'String',descr);
    set(mTextBox,'Units','inches');
    set(mTextBox,'Position',[3,5.5,3,2.5])

    subplot(3,2,3)
    plot(masterToneRaster{cellIndex(i,3)}(:,2),masterToneRaster{cellIndex(i,3)}(:,1),'b.')
    xlim(rasterWindow);
    ylim([1,length(toneTimes)]);
    set(gca, 'Units', 'inches');
    set(gca,'OuterPosition',[0.1 3 3 2.5]);
    title('Multiunit Raster Relative to Tone')

    subplot(3,2,4)
    plot(masterToneHist{cellIndex(i,3)}(:,2),masterToneHist{cellIndex(i,3)}(:,1),...
        masterToneHist{cellIndex(i,3)}(:,2),stePlotter(:,cellIndex(i,3),1),'b--',...
        masterToneHist{cellIndex(i,3)}(:,2),stePlotter(:,cellIndex(i,3),2),'b--')
    xlim(rasterWindow);
    set(gca, 'Units', 'inches');
    set(gca,'OuterPosition',[3.1 3 3 2.5]);
    title('Multiunit Histogram Relative to Tone')

    subplot(3,2,5)
    plot(masterToneRaster{cellIndex(i,2)}(:,2),masterToneRaster{cellIndex(i,2)}(:,1),'b.')
    xlim(rasterWindow);
    ylim([1,length(toneTimes)]);
    set(gca, 'Units', 'inches');
    set(gca,'OuterPosition',[0.1 0.2 3 2.5]);
    title('Unit Raster Relative to Tone')

    subplot(3,2,6)
    plot(masterToneHist{cellIndex(i,2)}(:,2),masterToneHist{cellIndex(i,2)}(:,1),...
        masterToneHist{cellIndex(i,2)}(:,2),stePlotter(:,cellIndex(i,2),1),'b--',...
        masterToneHist{cellIndex(i,2)}(:,2),stePlotter(:,cellIndex(i,2),2),'b--')
    xlim(rasterWindow);
    set(gca, 'Units', 'inches');
    set(gca,'OuterPosition',[3.1 0.2 3 2.5]);
    title('Unit Histogram Relative to Tone')
    
    



end



end