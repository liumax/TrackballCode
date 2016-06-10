%The purpose of this code is to extract rasters relative to a specific time
%point. This code also computs a histogram for all trials, and
%generates standard error lines for this.

function [structureTarget] = functionPairingRasterHistExtraction(i,clusters,...
    master,numDBs,numFreqs,structureTarget,truncatedNames,spikeTimes,...
    rasterWindow,histBin,histBinVector);
%holder for all full histograms
indivToneHist = cell(clusters,size(master,1)); 

%generate rasters and histograms
for j = 1:clusters
    rasterCounter = 1;
    for k = 1:size(master,1)
        spikeHolder = spikeTimes{j}...
            (spikeTimes{j}>master(k,1)+rasterWindow(1)...
            & spikeTimes{j}<master(k,1)+rasterWindow(2));
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
        counts = [];
        centers = [];
        %fills in large raster plot. holder updated position.
        toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,1) = k;
        toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,2) = spikeHolder;
        toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,3) = master(k,2); %stores frequency
        toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,4) = master(k,3); %stores amplitude
        toneRaster(rasterCounter:rasterCounter + length(spikeHolder)-1,5) = master(k,5); %stores freq/amp index
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

%the code below is for generating standard error lines for the plot of
%the overall histogram.
stdHolder = zeros(length(histBinVector),length(master(:,1)));
steHolder = zeros(length(histBinVector),clusters);

for j = 1:clusters
    for k = 1:length(master(:,1))
        stdHolder(:,k) = indivToneHist{j,k}(:,1);
    end
    steHolder(:,j) = std(stdHolder,0,2)/sqrt(length(master(:,1)));
end

%generates lines representing standard error.
stePlotter = zeros(length(histBinVector),clusters,2);
for j = 1:clusters
    stePlotter(:,j,1) = masterToneHist{j}(:,1)-steHolder(:,j);
    stePlotter(:,j,2) = masterToneHist{j}(:,1)+steHolder(:,j);
end

structureTarget.Rasters = masterToneRaster;
structureTarget.Histogram = masterToneHist;
structureTarget.StandardError = steHolder;
structureTarget.StandardErrorPlotting = stePlotter;

end






