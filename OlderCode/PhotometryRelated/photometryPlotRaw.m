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
    
    rawTimes = [0:1/s.Photo.RawRate:(1/s.Photo.RawRate)*(length(s.Photo.Raw)-1)];
    findZero = round(3/mean(diff(s.Photo.dFTime)));
    findEnd = round(4/mean(diff(s.Photo.dFTime)));
    %now lets make rasters of all the photometry data

    for i = 1:length(s.Photo.MBEDSig)
        targetTime = s.Photo.MBEDSig(i);
        rasterWindow = [-3,5];
        
        rasterSize = (rasterWindow(2)-rasterWindow(1))/mean(diff(s.Photo.dFTime));
        targetTime1 = find(s.Photo.dFTime - targetTime - rasterWindow(1)>0,1,'first');
        targetTime2 = round(targetTime1 + rasterSize);
        if targetTime2 < length(s.Photo.dFTime)
            rasterDF(:,i) = s.Photo.dFTrace(targetTime1:targetTime2);
            raster470(:,i) = s.Photo.x70(targetTime1:targetTime2);
            raster405(:,i) = s.Photo.x05(targetTime1:targetTime2);
        end
        
%         rasterSize2 = (rasterWindow(2)-rasterWindow(1))*s.Photo.RawRate;
%         targetTime1 = find(rawTimes - targetTime - rasterWindow(1)>0,1,'first');
%         targetTime2 = round(targetTime1 + rasterSize2);
%         
%         if targetTime2 < length(rawTimes)
%             rasterRaw(:,i) = downsample(s.Photo.Raw(3,targetTime1:targetTime2),100);
%         end
        
    end
    
    meanBigDF = mean(rasterDF(:,s.MBED.HiTrials)');
    meanSmallDF = mean(rasterDF(:,s.MBED.LowTrials)');
    
    meanBig470 = mean(raster470(:,s.MBED.HiTrials)');
    meanSmall470 = mean(raster470(:,s.MBED.LowTrials)');
    
    meanBig405 = mean(raster405(:,s.MBED.HiTrials)');
    meanSmall405 = mean(raster405(:,s.MBED.LowTrials)');
    
%     meanBigRaw = mean(rasterRaw(:,s.MBED.HiTrials)');
%     meanSmallRaw = mean(rasterRaw(:,s.MBED.LowTrials)');
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    subplot(2,2,1)
    hold on
    plot([findZero findZero],[min(min(raster405(:,s.MBED.HiTrials))) max(max(raster405(:,s.MBED.HiTrials)))],'r')
    plot([findEnd findEnd],[min(min(raster405(:,s.MBED.HiTrials))) max(max(raster405(:,s.MBED.HiTrials)))],'r')
    plot(raster405(:,s.MBED.HiTrials))
    title('405 Hi Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    subplot(2,2,3)
    hold on
    plot([findZero findZero],[min(min(raster405(:,s.MBED.LowTrials))) max(max(raster405(:,s.MBED.LowTrials)))],'r')
    plot([findEnd findEnd],[min(min(raster405(:,s.MBED.LowTrials))) max(max(raster405(:,s.MBED.LowTrials)))],'r')
    plot(raster405(:,s.MBED.LowTrials))
    title('405 Low Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    subplot(2,2,2)
    hold on
    plot([findZero findZero],[min(min(raster470(:,s.MBED.HiTrials))) max(max(raster470(:,s.MBED.HiTrials)))],'r')
    plot([findEnd findEnd],[min(min(raster470(:,s.MBED.HiTrials))) max(max(raster470(:,s.MBED.HiTrials)))],'r')
    plot(raster470(:,s.MBED.HiTrials))
    title('470 Hi Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    subplot(2,2,4)
    hold on
    plot([findZero findZero],[min(min(raster470(:,s.MBED.LowTrials))) max(max(raster470(:,s.MBED.LowTrials)))],'r')
    plot([findEnd findEnd],[min(min(raster470(:,s.MBED.LowTrials))) max(max(raster470(:,s.MBED.LowTrials)))],'r')
    plot(raster470(:,s.MBED.LowTrials))
    title('470 Low Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    hold off
    spikeGraphName = strcat('Individual Raw Traces',num2str(bigInd));
    savefig(hFig,spikeGraphName);
    
    
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    subplot(2,2,1)
    hold on
    plot([findZero findZero],[min(mean(raster405(:,s.MBED.HiTrials)')) max(mean(raster405(:,s.MBED.HiTrials)'))],'r')
    plot([findEnd findEnd],[min(mean(raster405(:,s.MBED.HiTrials)')) max(mean(raster405(:,s.MBED.HiTrials)'))],'r')
    plot(mean(raster405(:,s.MBED.HiTrials)'))
    title('Mean of 405 Hi Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    subplot(2,2,3)
    hold on
    plot([findZero findZero],[min(mean(raster405(:,s.MBED.LowTrials)')) max(mean(raster405(:,s.MBED.LowTrials)'))],'r')
    plot([findEnd findEnd],[min(mean(raster405(:,s.MBED.LowTrials)')) max(mean(raster405(:,s.MBED.LowTrials)'))],'r')
    plot(mean(raster405(:,s.MBED.LowTrials)'))
    title('Mean of 405 Low Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    subplot(2,2,2)
    hold on
    plot([findZero findZero],[min(mean(raster470(:,s.MBED.HiTrials)')) max(mean(raster470(:,s.MBED.HiTrials)'))],'r')
    plot([findEnd findEnd],[min(mean(raster470(:,s.MBED.HiTrials)')) max(mean(raster470(:,s.MBED.HiTrials)'))],'r')
    plot(mean(raster470(:,s.MBED.HiTrials)'))
    title('Mean of 470 Hi Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    subplot(2,2,4)
    hold on
    plot([findZero findZero],[min(mean(raster470(:,s.MBED.LowTrials)')) max(mean(raster470(:,s.MBED.LowTrials)'))],'r')
    plot([findEnd findEnd],[min(mean(raster470(:,s.MBED.LowTrials)')) max(mean(raster470(:,s.MBED.LowTrials)'))],'r')
    plot(mean(raster470(:,s.MBED.LowTrials)'))
    title('Mean of 470 Low Trials')
    xlabel('Samples')
    ylabel('Fluorescence (au)')
    xlim([1 length(raster405)])
    hold off
    spikeGraphName = strcat('Averaged Raw Traces',num2str(bigInd));
    savefig(hFig,spikeGraphName);
    
    
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    subplot(2,1,1)
    hold on
    plot((meanBigDF-min(meanBigDF))/(max(meanBigDF)-min(meanBigDF)),'k')
    plot((meanBig470-min(meanBig470))/(max(meanBig470)-min(meanBig470)),'g')
    plot((meanBig405-min(meanBig405))/(max(meanBig405)-min(meanBig405)),'b')
    xlim([1 length(raster405)])
    xlabel('Samples')
    title('Normalized Signals for Hi Trials, 405 (b) 470 (g) combined(k)')
    subplot(2,1,2)
    hold on
    plot((meanSmallDF-min(meanSmallDF))/(max(meanSmallDF)-min(meanSmallDF)),'k')
    plot((meanSmall470-min(meanSmall470))/(max(meanSmall470)-min(meanSmall470)),'g')
    plot((meanSmall405-min(meanSmall405))/(max(meanSmall405)-min(meanSmall405)),'b')
    xlim([1 length(raster405)])
    xlabel('Samples')
    title('Normalized Signals for Low Trials, 405 (b) 470 (g) combined(k)')
    hold off
    spikeGraphName = strcat('Normalized Traces',num2str(bigInd));
    savefig(hFig,spikeGraphName);
%     figure
%     subplot(2,1,1)
%     plot(meanBigRaw,'r')
%     subplot(2,1,2)
%     plot(meanSmallRaw,'r')
    
end