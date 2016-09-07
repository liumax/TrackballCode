%The purpose of this code is to extract rasters relative to a specific time
%point, and calculate average firing rate based on the period before the
%time point. This code also computs a histogram for all trials, and
%generates standard error lines for this.

%Main thing that needs explanation is how I'm calculating mean and standard
%deviation of baseline firing. Originally, I gathered every 5msec bin in
%the pre-tone period of 100msec, and then computed the mean and standard
%deviation based on this. The problem with this system was while the mean
%was correct, the standard deviation would be far too large to be possible
%(mean firing rate of 9 Hz, standard deviation of 42. or 0.75 Hz and 12
%STD). This is probably due to the fact that with so many 5msec bins, the
%standard deviation gets blown up massively. Instead, I am now averaging
%every trial to calculate a firing rate based on that trial. I then
%tabulate every single trial. The mean firing rate and standard deviation
%is computed off of these means. 

function [matclustStruct] = functionPulseBigRasterHistExtraction(i,clusterSizer,...
    master,baselineBins,matclustStruct,truncatedNames,...
    rasterWindow,histBin,histBinVector,laserSwitch);
%holder for all full histograms
indivToneHist = cell(clusterSizer,size(master,1)); 
%holds means of values from pre-tone period. Stores as array with columns = number of clusters
baselineHist = zeros(size(master,1),clusterSizer);

if laserSwitch == 1
    %% generate rasters and histograms for laser
    for j = 1:clusterSizer
        rasterCounter = 1;
        baselineCounter = 1;
        for k = 1:size(master,1)
            spikeHolder = matclustStruct.(truncatedNames{i}).SpikeTimes{j}...
                (matclustStruct.(truncatedNames{i}).SpikeTimes{j}>master(k,1)+rasterWindow(1)...
                & matclustStruct.(truncatedNames{i}).SpikeTimes{j}<master(k,1)+rasterWindow(2));
            spikeHolder = spikeHolder - master(k,1);
            %This will generate a histogram based on each individual trace,
            %which can be used to generate standard deviations.
            [counts,centers] = hist(spikeHolder,histBinVector);
            %below code necessary to prevent bugs with row vectors
            countSize = size(counts);
            centerSize = size(centers);
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            indivToneHist{j,k} = [counts'*(1/histBin),centers'];
            baselineHist(baselineCounter,j) = mean(counts(baselineBins(1):baselineBins(2)));
            baselineCounter = baselineCounter + 1;
            counts = [];
            centers = [];
            %fills in large raster plot. holder updated position.
            toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,1) = k;
            toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,2) = spikeHolder;
            toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,3) = master(k,2); %stores frequency
            if mod(k,2) == 0
                toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,4) = 1;
            else
                toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,4) = 0;
            end
            rasterCounter = rasterCounter + length(spikeHolder);
            toneRaster(toneRaster(:,1) == 0,:) = [];
            spikeHolder = [];
        end
        masterToneRaster{j} = toneRaster;
        toneRaster = zeros(100000,2);
        for k = 1:2
            [counts,centers] = hist(masterToneRaster{j}(masterToneRaster{j}(:,4) == (k-1),2),histBinVector);
            countSize = size(counts);
            centerSize = size(centers);
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            masterToneHist{j}{k} = [counts'*(1/histBin)/length(master(:,1)),centers'];
            counts = [];
            centers = [];
        end
    end
    %next I will calculate mean firing rate based on collected baseline histogram
    %information.

    preToneMean = mean(baselineHist)/histBin;
    preToneSTD = std(baselineHist)/histBin;%need hist bin division: std has NOT factored this in
    matclustStruct.(truncatedNames{i}).AverageFiringRate = preToneMean;
    matclustStruct.(truncatedNames{i}).AverageFiringSTD = preToneSTD;
    
    %need to adjust histograms to be z scored! Saves as new array to
    %preserve old data.

    masterZScore = masterToneHist;
    for j = 1:clusterSizer
        for k = 1:2
        masterZScore{j}{k}(:,1) = (masterToneHist{j}{k}(:,1)-preToneMean(j))/preToneSTD(j);
        end
    end
    
    %the code below is for generating standard error lines for the plot of
    %the overall histogram.
    stdHolder = zeros(length(histBinVector),length(master(:,1))/2,2);
    steHolder = zeros(length(histBinVector),clusterSizer,2);

    for j = 1:clusterSizer
        for k = 1:length(master(:,1))
            if mod(k,2) == 0
                stdHolder(:,k/2,2) = indivToneHist{j,k}(:,1);
            else
                stdHolder(:,(k+1)/2,1) = indivToneHist{j,k}(:,1);
            end
        end
        steHolder(:,j,1) = std(stdHolder(:,:,1),0,2)/sqrt(length(master(:,1))/2);
        steHolder(:,j,2) = std(stdHolder(:,:,2),0,2)/sqrt(length(master(:,1))/2);
    end
    
    %generates lines representing standard error.
    stePlotter = zeros(length(histBinVector),clusterSizer,2,2);
    for j = 1:clusterSizer
        stePlotter(:,j,1,1) = masterToneHist{j}{1}(:,1)-steHolder(:,j,1);
        stePlotter(:,j,2,1) = masterToneHist{j}{1}(:,1)+steHolder(:,j,1);
        stePlotter(:,j,1,2) = masterToneHist{j}{2}(:,1)-steHolder(:,j,2);
        stePlotter(:,j,2,2) = masterToneHist{j}{2}(:,1)+steHolder(:,j,2);
    end
else
    %% generate rasters and histograms for no laser
    for j = 1:clusterSizer
        rasterCounter = 1;
        baselineCounter = 1;
        for k = 1:size(master,1)
            spikeHolder = matclustStruct.(truncatedNames{i}).SpikeTimes{j}...
                (matclustStruct.(truncatedNames{i}).SpikeTimes{j}>master(k,1)+rasterWindow(1)...
                & matclustStruct.(truncatedNames{i}).SpikeTimes{j}<master(k,1)+rasterWindow(2));
            spikeHolder = spikeHolder - master(k,1);
            %This will generate a histogram based on each individual trace,
            %which can be used to generate standard deviations.
            [counts,centers] = hist(spikeHolder,histBinVector);
            %below code necessary to prevent bugs with row vectors
            countSize = size(counts);
            centerSize = size(centers);
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            indivToneHist{j,k} = [counts'*(1/histBin),centers'];
            baselineHist(baselineCounter,j) = mean(counts(baselineBins(1):baselineBins(2)));
            baselineCounter = baselineCounter + 1;
            counts = [];
            centers = [];
            %fills in large raster plot. holder updated position.
            toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,1) = k;
            toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,2) = spikeHolder;
            toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,3) = master(k,2); %stores frequency
            rasterCounter = rasterCounter + length(spikeHolder);
            toneRaster(toneRaster(:,1) == 0,:) = [];
            spikeHolder = [];
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
        masterToneHist{j} = [counts'*(1/histBin)/length(master(:,1)),centers'];
        counts = [];
        centers = [];
    end
    %next I will calculate mean firing rate based on collected baseline histogram
    %information.

    preToneMean = mean(baselineHist)/histBin;
    preToneSTD = std(baselineHist)/histBin;%need hist bin division: std has NOT factored this in
    matclustStruct.(truncatedNames{i}).AverageFiringRate = preToneMean;
    matclustStruct.(truncatedNames{i}).AverageFiringSTD = preToneSTD;

    %need to adjust histograms to be z scored! Saves as new array to
    %preserve old data.

    masterZScore = masterToneHist;
    for j = 1:clusterSizer
        masterZScore{j}(:,1) = (masterToneHist{j}(:,1)-...
            matclustStruct.(truncatedNames{i}).AverageFiringRate(j))...
            /matclustStruct.(truncatedNames{i}).AverageFiringSTD(j);
    end

    % masterZScore = masterToneHist;
    % for j = 1:clusterSizer
    %     masterZScore{j}(:,1) = zscore(masterToneHist{j}(:,1));
    % end

    %the code below is for generating standard error lines for the plot of
    %the overall histogram.
    stdHolder = zeros(length(histBinVector),length(master(:,1)));
    steHolder = zeros(length(histBinVector),clusterSizer);

    for j = 1:clusterSizer
        for k = 1:length(master(:,1))
            stdHolder(:,k) = indivToneHist{j,k}(:,1);
        end
        steHolder(:,j) = std(stdHolder,0,2)/sqrt(length(master(:,1)));
    end

    %generates lines representing standard error.
    stePlotter = zeros(length(histBinVector),clusterSizer,2);
    for j = 1:clusterSizer
        stePlotter(:,j,1) = masterToneHist{j}(:,1)-steHolder(:,j);
        stePlotter(:,j,2) = masterToneHist{j}(:,1)+steHolder(:,j);
    end
    
end

matclustStruct.(truncatedNames{i}).BigRasters = masterToneRaster;
matclustStruct.(truncatedNames{i}).BigHistogram = masterToneHist;
matclustStruct.(truncatedNames{i}).PreToneFiringRate = preToneMean;
matclustStruct.(truncatedNames{i}).PreToneFiringSTD = preToneSTD;
matclustStruct.(truncatedNames{i}).BigZScore = masterZScore;
matclustStruct.(truncatedNames{i}).BigStandardError = steHolder;
matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting = stePlotter;

end






