%This is dirty code to re-run all the plotting I did with previous things
%using just the raw data from 470 and 405 feeds. 

testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
targetNames = testNames(findString);

numFiles = length(targetNames);

for bigInd = 1:numFiles
    %pull file and load
    fileName = targetNames{bigInd};
    load(fileName)
    
    %find time differences between data points
    timeDiff = mean(diff(s.Photo.dFTime));
    timeVector = [-2:timeDiff:3];
    findZero = round(2/timeDiff);
    findEnd = round(2.1/timeDiff);
    %now lets make rasters of all the photometry data

    for i = 1:length(s.Photo.AlignTimes)
        targetTime = s.Photo.AlignTimes(i);
        rasterWindow = [-2,3];
        
        rasterSize = (rasterWindow(2)-rasterWindow(1))/mean(diff(s.Photo.dFTime));
        targetTime1 = find(s.Photo.dFTime - targetTime - rasterWindow(1)>0,1,'first');
        targetTime2 = round(targetTime1 + rasterSize);
        if targetTime2 < length(s.Photo.dFTime)
            rasterDF(:,i) = s.Photo.dFTrace(targetTime1:targetTime2);
            raster470(:,i) = s.Photo.x70(targetTime1:targetTime2);
            raster405(:,i) = s.Photo.x05(targetTime1:targetTime2);
        end
        
    end
    
    overMeanDF = mean(rasterDF');
    overMean470 = mean(raster470');
    overMean405 = mean(raster405');
    
    %now pull sound data to isolate by amplitude
    
    amps = s.SoundData.Amplitudes;
    uniqueAmps = unique(amps);
    ampInds = zeros(length(amps)/3,3);
    
    freqs = s.SoundData.Frequencies;
    uniqueFreqs = unique(freqs);
    indStore = [];
    for i = 1:length(uniqueAmps)
        for j = 1:length(uniqueFreqs)
            indStore(i,j,:) = find(amps == uniqueAmps(i) & freqs == uniqueFreqs(j));
        end
    end
    
    maxAmps = find(amps == uniqueAmps(end));
    minAmps = find(amps == uniqueAmps(1));
    
    maxAmpDF = mean(rasterDF(:,maxAmps)');
    maxAmp470 = mean(raster470(:,maxAmps)');
    maxAmp405 = mean(raster405(:,maxAmps)');
    
    minAmpDF = mean(rasterDF(:,minAmps)');
    minAmp470 = mean(raster470(:,minAmps)');
    minAmp405 = mean(raster405(:,minAmps)');
    
    %make heatmaps by frequency and DB
    for i = 1:length(uniqueAmps)
        for j = 1:length(uniqueFreqs)
            meanStoreDF(i,j,:) = mean(rasterDF(:,indStore(i,j,:))');
            meanStore470(i,j,:) = mean(raster470(:,indStore(i,j,:))');
            meanStore405(i,j,:) = mean(raster405(:,indStore(i,j,:))');
        end
    end
    
    limsDF = [min(min(min(meanStoreDF))) max(max(max(meanStoreDF)))];
    lims470 = [min(min(min(meanStore470))) max(max(max(meanStore470)))];
    lims405 = [min(min(min(meanStore405))) max(max(max(meanStore405)))];
    
%     meanBigRaw = mean(rasterRaw(:,s.MBED.HiTrials)');
%     meanSmallRaw = mean(rasterRaw(:,s.MBED.LowTrials)');


    %create a set of points for the raster axis
    rasterAxis = zeros(rasterWindow(2)-rasterWindow(1)+1,2);
    rasterAxis(:,1) = [rasterWindow(1):rasterWindow(2)];
    %find a second in photometry sample time
    rasterAxis(:,2) = [0:length(minAmpDF)/(rasterWindow(2) - rasterWindow(1)):length(minAmpDF)];
    rasterAxis(1,2) = 1;
    
    octaveRange = zeros(totalOctaves + 1,2);
    octaveRange(1,1) = uniqueFreqs(1);
    for ind = 1:totalOctaves
        octaveRange (ind+1,1) = octaveRange(ind,1)*2;
    end
    %next, I find the positions from uniqueFreqs that match octaveRange
    for ind = 1:size(octaveRange,1);
        octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
    end


    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    subplot(3,4,1)
    hold on
    plot(timeVector,maxAmp405)
    plot(timeVector,minAmp405,'r')
    plot([0 0],[min(maxAmp405) max(maxAmp405)],'k')
    plot([s.SoundData.ToneDuration s.SoundData.ToneDuration],[min(maxAmp405) max(maxAmp405)],'k')
    xlim(rasterWindow)
    title('405 max(b) and min(r) amp')
    
    subplot(3,4,2)
    imagesc(squeeze(meanStore405(1,:,:)),lims405)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title(fileName)
    
    subplot(3,4,3)
    imagesc(squeeze(meanStore405(2,:,:)),lims405)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    subplot(3,4,4)
    imagesc(squeeze(meanStore405(3,:,:)),lims405)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    subplot(3,4,5)
    plot(timeVector,maxAmp470)
    hold on
    plot(timeVector,minAmp470,'r')
    plot([0 0],[min(maxAmp470) max(maxAmp470)],'k')
    plot([s.SoundData.ToneDuration s.SoundData.ToneDuration],[min(maxAmp470) max(maxAmp470)],'k')
    xlim(rasterWindow)
    title('470 max(b) and min(r) amp')
    
    subplot(3,4,6)
    imagesc(squeeze(meanStore470(1,:,:)),lims470)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    subplot(3,4,7)
    imagesc(squeeze(meanStore470(2,:,:)),lims470)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    subplot(3,4,8)
    imagesc(squeeze(meanStore470(3,:,:)),lims470)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    subplot(3,4,9)
    plot(timeVector,maxAmpDF)
    hold on
    plot(timeVector,minAmpDF,'r')
    plot([0 0],[min(maxAmpDF) max(maxAmpDF)],'k')
    plot([s.SoundData.ToneDuration s.SoundData.ToneDuration],[min(maxAmpDF) max(maxAmpDF)],'k')
    xlim(rasterWindow)
    title('DF max(b) and min(r) amp')
    
    subplot(3,4,10)
    imagesc(squeeze(meanStoreDF(1,:,:)),limsDF)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    subplot(3,4,11)
    imagesc(squeeze(meanStoreDF(2,:,:)),limsDF)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    subplot(3,4,12)
    imagesc(squeeze(meanStoreDF(3,:,:)),limsDF)
    colormap parula
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    
    
    

    spikeGraphName = strcat('Individual Raw Traces',num2str(bigInd));
    savefig(hFig,spikeGraphName);
    
    
end