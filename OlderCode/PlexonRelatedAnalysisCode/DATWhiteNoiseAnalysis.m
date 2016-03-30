

fileName = '5ms pulse 20 hz 20 pulses 20isi'

%Reads in NEX File
[nexFile] = readNexFile(strcat(char(fileName),'.nex'));

%Variables I may want to adjust
rasterWindow = [-1,2];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
histBin = 0.1; %bin size in seconds
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
%determines how many repeated events there were
realEvents = find(eventTstampSize>1);
realEventNum = length(realEvents);
%% 

if realEventNum ==1 %if only a single event!
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
    
    plotSizer = length(masterToneHist);

    figure
    for i = 1:plotSizer
        subplot(2,plotSizer,i)
        plot(masterToneRaster{i}(:,2),masterToneRaster{i}(:,1),'b.')
    %     title(['Cell #',num2str(i),' Raster'])
    %     xlabel('Seconds')
        xlim(rasterWindow);
        subplot(2,plotSizer,i+plotSizer)
        plot(masterToneHist{i}(:,2),masterToneHist{i}(:,1),...
            masterToneHist{i}(:,2),stePlotter(:,i,1),'b--',...
            masterToneHist{i}(:,2),stePlotter(:,i,2),'b--')
    %     title(['Cell #',num2str(i),' Histogram with Binsize ',num2str(histBin)])
    %     xlabel('Seconds')
    %     ylabel('Average Firing Rate (Hz)')
        xlim(rasterWindow);
    end
    %% 
elseif realEventNum == 2 
    %This code relies on assumption that stim trials are fewer than
    %non-stim trials.
    bigEvents = find(eventTstampSize == max(eventTstampSize));
    smallEvents = find(realEvents ~= bigEvents);
    toneTimes = nexFile.events{bigEvents}.timestamps;
    toneStimTimes= nexFile.events{smallEvents}.timestamps;
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
    rasterHolder = 1;
    toneRaster = zeros(100000,2);
    masterToneRaster = cell(cellSize,1);
    masterToneHist = cell(cellSize,1);
    indivToneHist = cell(cellSize,length(toneTimes));
    for j = 1:cellSize
        for i = 1:length(toneTimes)
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
    
    
    rasterHolder = 1;
    toneStimRaster = zeros(100000,2);
    masterToneStimRaster = cell(cellSize,1);
    masterToneStimHist = cell(cellSize,1);
    indivToneStimHist = cell(cellSize,length(toneStimTimes));
    for j = 1:cellSize
        for i = 1:length(toneStimTimes)
            toneHolder = cellData{j}(cellData{j}>toneStimTimes(i)+rasterWindow(1) & cellData{j} < toneStimTimes(i) + rasterWindow(2));
            toneHolder = toneHolder - toneStimTimes(i);
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
            indivToneStimHist{j,i} = [counts'*(1/histBin),centers'];
            counts = [];
            centers = [];
            toneStimRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
            toneStimRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
            rasterHolder = rasterHolder + length(toneHolder);
            toneStimRaster(toneStimRaster(:,1) == 0,:) = [];
            toneHolder = [];
        end
        masterToneStimRaster{j} = toneStimRaster;
        toneStimRaster = zeros(100000,2);
        [counts,centers] = hist(masterToneStimRaster{j}(:,2),histBinVector);
        countSize = size(counts);
        centerSize = size(centers);
        if countSize(1)>countSize(2)
            counts = counts';
        end
        if centerSize(1)>centerSize(2)
            centers = centers';
        end
        masterToneStimHist{j} = [counts'*(1/histBin)/length(toneStimTimes),centers'];
        counts = [];
        centers = [];
    end
    
    stdStimHolder = zeros(length(histBinVector),length(toneStimTimes));
    steStimHolder = zeros(length(histBinVector),cellSize);

    for i = 1:cellSize
        for j = 1:length(toneStimTimes)
            stdStimHolder(:,j) = indivToneStimHist{i,j}(:,1);
        end
        steStimHolder(:,i) = std(stdStimHolder,0,2)/sqrt(length(toneStimTimes));
    end
    
    steStimPlotter = zeros(length(histBinVector),cellSize,2);
    for i = 1:cellSize
        steStimPlotter(:,i,1) = masterToneStimHist{i,1}(:,1)-steStimHolder(:,i);
        steStimPlotter(:,i,2) = masterToneStimHist{i,1}(:,1)+steStimHolder(:,i);
    end
    
    %Adjusts rasters so that all of them are aligned to the appropriate
    %order. This will allow things to be properly interleaved when
    %displayed in the raster.
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

    figure
    for i = 1:plotSizer
        subplot(2,plotSizer,i)
        plot(masterToneRaster{i}(:,2),masterToneRaster{i}(:,1),'b.',...
            masterToneStimRaster{i}(:,2),masterToneStimRaster{i}(:,1),'r.')
    %     title(['Cell #',num2str(i),' Raster'])
    %     xlabel('Seconds')
        subplot(2,plotSizer,i+plotSizer)
        plot(masterToneHist{i}(:,2),masterToneHist{i}(:,1),...
            masterToneHist{i}(:,2),stePlotter(:,i,1),'b--',...
             masterToneHist{i}(:,2),stePlotter(:,i,2),'b--',...
            masterToneStimHist{i}(:,2),masterToneStimHist{i}(:,1),'r',...
             masterToneStimHist{i}(:,2),steStimPlotter(:,i,1),'r--',...
             masterToneStimHist{i}(:,2),steStimPlotter(:,i,2),'r--')
    %     title(['Cell #',num2str(i),' Histogram with Binsize ',num2str(histBin)])
    %     xlabel('Seconds')
    %     ylabel('Average Firing Rate (Hz)')
        xlim(rasterWindow);
    end
end



