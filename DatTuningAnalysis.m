%Variables I may want to adjust
rasterWindow = [-1,2];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
histBin = 0.05; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)-histBin/2:histBin:rasterWindow(2)-histBin/2];
eventNum = 1; %This is the number for the event number you wish to use. 
%For my current experiments, 1 is for pulseblaster/NIDAQ, 2 is for mBED
colorSteps = 5; %This determines how many steps you take in changing colors. 
%5 means that for RGB, each value will change 0.2 every step
soundDur = 500; %sound duration, in MS


%insert matlab file name here
matName = 'ML150730B_150805_AP14DV3468_soundSet'
%open matlab file
load(matName);
toneRecord(toneRecord(:,2)==0,:) = [];

fileName = 'ML150730B_AP14_DV3468_tuning'

%Reads in NEX File
[nexFile] = readNexFile(strcat(char(fileName),'.nex'));

%picks out correct events to set tone times
eventSize = length(nexFile.events);
eventNames = cell(eventSize,1);
toneTimes= nexFile.events{eventNum}.timestamps;
toneTimes(1) = [];

%readjusts tone times for toneRecord
toneRecord(:,1)= toneTimes;

%Preps cell array for firing data.
cellSize = length(nexFile.neurons);
cellData = cell(cellSize,1);

for i = 1:cellSize
    cellData{i} = nexFile.neurons{i}.timestamps;
end


%Generates raster plot for tone times as well as histogram. Also retains
%frequency information
rasterHolder = 1;
toneRaster = zeros(100000,3);
masterToneRaster = cell(cellSize,1);
masterToneHist = cell(cellSize,1);
for j = 1:cellSize
    for i = 1:length(toneTimes)
        toneHolder = cellData{j}(cellData{j}>toneTimes(i)+rasterWindow(1) & cellData{j} < toneTimes(i) + rasterWindow(2));
        toneHolder = toneHolder - toneTimes(i);
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,3) = toneRecord(i,2);
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


% tuningRaster = masterToneRaster;
% for j = 1:length(tuningRaster)
%     for i = 1:length(toneRecord)
%         tuningRaster{j}(tuningRaster{j}(:,1)==i,1) = toneRecord(i,2);
%     end
% end

%processes tone value information
toneValues = union(toneRecord(:,2),toneRecord(:,2));
toneReps = zeros(length(toneValues),2);
for i = 1:length(toneValues)
    toneReps(i,1) = toneValues(i);
    toneReps(i,2) = length(find(toneRecord(:,2) == toneValues(i)));
end
tuningSpecRaster = cell(length(toneValues),length(masterToneRaster));
tuningSpecHist = cell(length(toneValues),length(masterToneRaster));

%generates cell array with all data points associated with a specific
%frequency. columns are cells, rows are frequencies
for j = 1:length(masterToneRaster)
    for i = 1:length(toneValues)
        tuningSpecRaster{i,j}= masterToneRaster{j}(masterToneRaster{j}(:,3)==toneValues(i),2);
        [counts,centers] = hist(tuningSpecRaster{i,j},histBinVector);
        countSize = size(counts);%this is here because matlab is dumb and sometimes it comes as columns
        centerSize = size(centers); %sometimes as rows
        if countSize(1)>countSize(2) && centerSize(1) > centerSize(2)
            tuningSpecHist{i,j} = [counts,centers];
        elseif countSize(1)<countSize(2) && centerSize(1) < centerSize(2)
            tuningSpecHist{i,j} = [counts',centers'];
        elseif countSize(1)>countSize(2) && centerSize(1) < centerSize(2)
            tuningSpecHist{i,j} = [counts,centers'];
        elseif countSize(1)<countSize(2) && centerSize(1) > centerSize(2)
            tuningSpecHist{i,j} = [counts',centers];
        end
        
    end
end

%generate array for choice of color
colorArray = cell(length(toneValues),1);

for i = 1:length(colorArray)
    if i <= colorSteps
        colorArray{i} = [0 0 i/colorSteps];
    elseif i > colorSteps & i <= 2*colorSteps
        colorArray{i} = [(i-colorSteps)/colorSteps 0  1-((i-colorSteps)/colorSteps)];
    elseif i > 2*colorSteps & i <= 3*colorSteps
        colorArray{i} = [1-(i-2*colorSteps)/colorSteps (i-2*colorSteps)/colorSteps 0];
    elseif i > length(colorArray)
        disp('FAILURE OF COLOR SPECTRUM')
    end
end

[counts,centers] = hist(masterToneRaster{j}(:,2),histBinVector);

plotSizer = length(masterToneHist);

figure
for i = 1:plotSizer
    subplot(3,plotSizer,i)
    plot(masterToneRaster{i}(:,2),masterToneRaster{i}(:,1),'k.')
    title(['Cell #',num2str(i),' Raster'])
    xlabel('Seconds')
    xlim(rasterWindow)
    subplot(3,plotSizer,i+plotSizer)
    plot(masterToneHist{i}(:,2),masterToneHist{i}(:,1)*10/length(toneTimes),'k')
    title(['Cell #',num2str(i),' Hist Binsize ',num2str(histBin)])
    xlabel('Seconds')
    ylabel('Av. Firing Rate (Hz)')
    xlim(rasterWindow);
    subplot(3,plotSizer,i+2*(plotSizer))
    hold on
    for j = 1:length(toneValues)
        plot(tuningSpecHist{j,i}(:,2),tuningSpecHist{j,i}(:,1)*10/toneReps(i,2),...
            'Color',colorArray{j})
    end
    hold off
    xlim(rasterWindow);
    xlabel('Seconds')
    ylabel('Av. Firing Rate (Hz)')
    title(['Cell #',num2str(i),' Histogram by Tone Freq'])
end














