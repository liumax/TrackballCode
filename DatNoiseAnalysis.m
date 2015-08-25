%This is code meant for the analysis of just noise presentation, without
%accompanying dopamine stimulation or tuning curves

fileName = 'ML150730_AP17_DV2070_3to6secISIrand_NoDAT'

%Reads in NEX File
[nexFile] = readNexFile(strcat(char(fileName),'.nex'));

%Variables I may want to adjust
rasterWindow = [-1,2];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
histBin = 0.05; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)-histBin/2:histBin:rasterWindow(2)-histBin/2];

%Currently not in use, could be used to automate event name detection
eventSize = length(nexFile.events);
eventNames = cell(eventSize,1);

%Fills in Tone and Tone+Stim times, also figures out how Tone+Stim
%integrates with regular Tone delivery.
toneTimes= nexFile.events{2}.timestamps;
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

plotSizer = length(masterToneHist);

figure
for i = 1:plotSizer
    subplot(2,plotSizer,i)
    plot(masterToneRaster{i}(:,2),masterToneRaster{i}(:,1),'b.')
    title(['Cell #',num2str(i),' Raster'])
    xlabel('Seconds')
    xlim(rasterWindow)
    subplot(2,plotSizer,i+plotSizer)
    plot(masterToneHist{i}(:,2),masterToneHist{i}(:,1)*(1/histBin)/length(toneTimes))
    title(['Cell #',num2str(i),' Histogram with Binsize ',num2str(histBin)])
    xlabel('Seconds')
    ylabel('Average Firing Rate (Hz)')
    xlim(rasterWindow);
end


