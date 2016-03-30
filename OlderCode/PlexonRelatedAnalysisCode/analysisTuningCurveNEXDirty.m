%this is to be analysis code for tuning curve data!
rasterWindow = [0,0.05];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];

eventNum = 2; %This is the number for the event number you wish to use. 
%For my current experiments, 1 is for pulseblaster/NIDAQ, 2 is for mBED

%%
%NEED FILE INPUT FOR MATLAB FILE (TUNING CODE)
matName = '151209_ML151202A_R12_2500_tuning';
load(matName);

freqs = soundData.Frequencies;
uniqueFreqs = unique(freqs);
dBs = soundData.dBs;
uniqueDBs = unique(dBs);
Amplitudes = soundData.Amplitudes;

[sortedFreq,freqIndex] = sort(freqs);

%%
%NEED INPUT FROM PLEXON
nexName = '151209_ML1501202A_R12_2500_tuning';
[nexFile] = readNexFile(strcat(char(nexName),'.nex'));


%%


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

%%
%picks out correct events to set tone times
eventSize = length(nexFile.events);
toneTimes= nexFile.events{eventNum}.timestamps;

eventNames = cell(eventSize,1);
eventTstamps = cell(eventSize,1);
eventTstampSize = zeros(eventSize,1);

for i= 1:eventSize
    eventNames{i} = nexFile.events{i,1}.name;
    eventTstamps{i} = nexFile.events{i,1}.timestamps;
    eventTstampSize(i) = length(eventTstamps{i});
end

%Preps cell array for firing data.
cellSize = length(nexFile.neurons);
cellData = cell(cellSize,1);

for i = 1:cellSize
    cellData{i} = nexFile.neurons{i}.timestamps;
end

%%
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

%%
%This is to compute average firing rate of every 'true' unit.
averageFiringRate = zeros(trueUnits,1);
for i = 1:trueUnits
    averageFiringRate(i) = length(nexFile.neurons{unitCells(i),1}.timestamps)/totalTime;
end

%%
%Generates raster plot for tone times as well as histogram. Also retains
%frequency information
rasterHolder = 1;
toneRaster = zeros(100000,4);
masterToneRaster = cell(cellSize,1);

for j = 1:cellSize
    for i = 1:length(toneTimes)
        toneHolder = cellData{j}(cellData{j}>toneTimes(i)+rasterWindow(1) & cellData{j} < toneTimes(i) + rasterWindow(2));
        toneHolder = toneHolder - toneTimes(i);
        %fills in large raster plot. holder updated position.
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,3) = freqs(i,1);
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,4) = dBs(i,1);
        rasterHolder = rasterHolder + length(toneHolder);
        toneRaster(toneRaster(:,1) == 0,:) = [];
        toneHolder = [];
    end
    masterToneRaster{j} = toneRaster;
    toneRaster = zeros(100000,2);
end

%this will convert the overall data into data specific for specific
%frequencies and amplitudes. Each cell represents information from a single
%channel. The rows represent amplitude, the columns represent frequency.
%Counts all spikes within the tone period (I know this is shitty)
processRaster = cell(cellSize,1);

for i = 1:cellSize
    processRaster{i} = zeros(length(uniqueDBs),length(uniqueFreqs));
    for j = 1:length(uniqueFreqs)
        for k = 1:length(uniqueDBs)
            processRaster{i}(k,j) = length(find(masterToneRaster{i}(:,3) == uniqueFreqs(j) & masterToneRaster{i}(:,4) == uniqueDBs(k)));
        end
    end
end

% this will code the image using matlabs built in color system. should make my own eventually.
% image(processRaster{1})
% colorbar

%%
%actual graphing code
set(0, 'DefaulttextInterpreter', 'none')
for i = 1:length(cellIndex)
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
        nexName;
        'Unit:';
        unitNames{cellIndex(i,2)}
        'Average Firing Rate (Hz):';
        averageFiringRate(i)};
    set(mTextBox,'String',descr);
    set(mTextBox,'Units','inches');
    set(mTextBox,'Position',[3,5.5,3,2.5])
    
    subplot(2,2,3)
    image(processRaster{cellIndex(i,3)})
    title('Tuning Curve of Multiunit')
    set(gca,'YTickLabel',uniqueDBs,...
        'XTickLabel',uniqueFreqs)
    
    subplot(2,2,4)
    image(processRaster{cellIndex(i,2)})
    set(gca,'YTickLabel',uniqueDBs,...
        'XTickLabel',uniqueFreqs)
    colorbar
    title('Tuning Curve of Single Unit')
end







