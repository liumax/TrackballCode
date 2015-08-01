

fileName = 'Tone and DA Trials DV34'

%Reads in NEX File
[nexFile] = readNexFile(strcat(char(fileName),'.nex'));

%Variables I may want to adjust
rasterWindow = [-1,2];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
histBin = 0.1; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)-histBin/2:histBin:rasterWindow(2)-histBin/2];

%Currently not in use, could be used to automate event name detection
eventSize = length(nexFile.events);
eventNames = cell(eventSize,1);

%Fills in Tone and Tone+Stim times, also figures out how Tone+Stim
%integrates with regular Tone delivery.
toneTimes= nexFile.events{1}.timestamps;
toneStimTimes= nexFile.events{2}.timestamps;
fullTimes = union(toneTimes,toneStimTimes);
[y,ia,ib]=intersect(fullTimes,toneStimTimes);
toneStimTimesIndex = ia;
[y,ia,ib]=intersect(fullTimes,toneTimes);
toneTimesIndex = ia;

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
for j = 1:cellSize
    for i = 1:length(toneTimes)
        toneHolder = cellData{j}(cellData{j}>toneTimes(i)+rasterWindow(1) & cellData{j} < toneTimes(i) + rasterWindow(2));
        toneHolder = toneHolder - toneTimes(i);
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
        rasterHolder = rasterHolder + length(toneHolder);
        toneRaster(toneRaster(:,1) == 0,:) = [];
        toneHolder = [];
    end
    masterToneRaster{j} = toneRaster;
    toneRaster = zeros(100000,2);
    [counts,centers] = hist(masterToneRaster{j}(:,2),histBinVector);
    masterToneHist{j} = [counts',centers'];
    counts = [];
    centers = [];
end

%Generates raster plot and histogram for tone+stim times
rasterHolder = 1;
toneStimRaster = zeros(100000,2);
masterToneStimRaster = cell(cellSize,1);
masterToneStimHist = cell(cellSize,1);
for j = 1:cellSize
    for i = 1:length(toneStimTimes)
        toneHolder = cellData{j}(cellData{j}>toneStimTimes(i)+rasterWindow(1) & cellData{j} < toneStimTimes(i) + rasterWindow(2));
        toneHolder = toneHolder - toneStimTimes(i);
        toneStimRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
        toneStimRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
        rasterHolder = rasterHolder + length(toneHolder);
        toneStimRaster(toneStimRaster(:,1) == 0,:) = [];
        toneHolder = [];
    end
    masterToneStimRaster{j} = toneStimRaster;
    toneStimRaster = zeros(100000,2);
    [counts,centers] = hist(masterToneStimRaster{j}(:,2),histBinVector);
    masterToneStimHist{j} = [counts',centers'];
    counts = [];
    centers = [];
end

%Adjusts rasters so that all of them are aligned by actual time of delivery
for j = 1:length(masterToneRaster)
    for i = length(toneTimesIndex):-1:1
        x=find(masterToneRaster{j}(:,1) == i);
        masterToneRaster{j}(x,1) = toneTimesIndex(i);
    end
end

for j = 1:length(masterToneStimRaster)
    for i = length(toneStimTimesIndex):-1:1
        x=find(masterToneStimRaster{j}(:,1) == i);
        masterToneStimRaster{j}(x,1) = toneStimTimesIndex(i);
    end
end

plotSizer = length(masterToneHist);

figure('Name','Blue is Tone Only, Red is Tone + Stim','NumberTitle','off')
for i = 1:plotSizer
    subplot(2,plotSizer,i)
    plot(masterToneRaster{i}(:,2),masterToneRaster{i}(:,1),'b.',...
        masterToneStimRaster{i}(:,2),masterToneStimRaster{i}(:,1),'r.')
    title(['Cell #',num2str(i),' Raster'])
    xlabel('Seconds')
    subplot(2,plotSizer,i+plotSizer)
    plot(masterToneHist{i}(:,2),masterToneHist{i}(:,1)*10/length(toneTimes),...
        masterToneStimHist{i}(:,2),masterToneStimHist{i}(:,1)*10/length(toneStimTimes),'r')
    title(['Cell #',num2str(i),' Histogram with Binsize ',num2str(histBin)])
    xlabel('Seconds')
    ylabel('Average Firing Rate (Hz)')
    xlim(rasterWindow);
end


