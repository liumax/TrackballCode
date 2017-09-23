%This is meant to be wrapper function for analysis of tuning curve
%photometry data.


targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);

for bigInd = 1:numFiles
    
    %load the target file
    load(targetFiles{bigInd})
    
    fullMax = max(s.Photo.Photo.x70dF);
    
    
    photoTimeStep = mean(diff(s.Photo.Photo.x70dFTime));
    rasterPhotWindow = round([-3,5]/photoTimeStep);
    for ind = 1:length(s.Photo.MBEDSig)
        alignTime = s.Photo.MBEDSig(ind);
        %find the time in the photometry trace
        photoPoint = find(s.Photo.Photo.x70dFTime - alignTime > 0,1,'first');
        if photoPoint + rasterPhotWindow(2) < length(s.Photo.Photo.x70dF)
            photoRaster(:,ind) = s.Photo.Photo.x70dF(photoPoint + rasterPhotWindow(1):photoPoint + rasterPhotWindow(2));
        else
            disp('Tone Rasters: Reached End of Photometry Trace')
            disp(ind)
            break
        end
    end
    rasterVect = [-3:photoTimeStep:5];
    zeroPoint = find(rasterVect > 0,1,'first');
    endPoint = find(rasterVect > 1,1,'first');
    %baseline period will simply be average of the five bins before zero
    zeroBase = mean(photoRaster(zeroPoint - 4:zeroPoint,:));
    endBase =  mean(photoRaster(endPoint - 4:endPoint,:));
    %find peak values for the twentyfive bins following the stimulus onset
    %(~500ms)
    zeroPeak = max(photoRaster(zeroPoint:zeroPoint+ 25,:));
    endPeak = max(photoRaster(zeroPoint:zeroPoint+ 25,:));
    
    zeroPeakVal = zeroPeak - zeroBase;
    endPeakVal = endPeak-endBase;
    

    numFreqs = length(unique(s.SoundData.Frequencies));
    numDBs = length(unique(s.SoundData.dBs));
    uniqueFreqs = unique(s.SoundData.Frequencies);
    uniqueDBs = unique(s.SoundData.dBs);



    trialMatrix = s.SoundData.TrialMatrix;
    for ind = 1:numFreqs
        for j = 1:numDBs
            %find all trials of the particular setting
            targetFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == uniqueDBs(j));
            %pull and average these traces
            tempHolder = photoRaster(:,targetFinder);
            photoRasterStore{ind,j} = tempHolder;
            peakValStore{ind,j} = (zeroPeakVal(targetFinder));
            endPeakValStore{ind,j} = (endPeakVal(targetFinder));
            peakValAV(ind,j) = mean(zeroPeakVal(targetFinder));
            endPeakValAV(ind,j) = mean(endPeakVal(targetFinder));
            tempHolder = mean(tempHolder');
            photoAverages(:,ind,j) = tempHolder;
        end
    end

    subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    subplot(2,2,1)
    imagesc(squeeze(photoAverages(:,:,3)'))
    colormap('parula')
%     colormap('jet')
    colorbar
    title('Responses over time, frequency increasing as you go down')
    
    subplot(2,2,3)
    hist(s.Photo.Peaks(:,1),100)
    title(targetFiles{bigInd})
    
    
    c1 = min(min([peakValAV,endPeakValAV]));
    c2 = max(max([peakValAV,endPeakValAV]));
    
    subplot(2,2,2)
    imagesc(peakValAV,[c1,c2])
    colormap('parula')
    colorbar
    title('Averaged Onset Responses')
    
    subplot(2,2,4)
    imagesc(endPeakValAV,[c1,c2])
    colormap('parula')
    colorbar
    title('Averaged Offset Responses')
    
    spikeGraphName = strcat((targetFiles{bigInd}(1:end-4)),'Figure');

    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
    fillName = (targetFiles{bigInd}(1:end-4));
    medVal = median(s.Photo.Peaks(:,1));
    
    fullStore.(fillName).peakValAV = peakValAV;
    fullStore.(fillName).endPeakValAV = endPeakValAV;
    fullStore.(fillName).peakValAVCorr = peakValAV/medVal;
    fullStore.(fillName).endPeakValAVCorr = endPeakValAV/medVal;
    fullStore.(fillName).peakValStore = peakValStore;
    fullStore.(fillName).endPeakValStore = endPeakValStore;
    fullStore.(fillName).photoAverages = photoAverages;
    fullStore.(fillName).peakVals = s.Photo.Peaks(:,1);
    fullStore.(fillName).medVal = medVal;
    fullStore.Names{bigInd} = fillName;
    
    
end



%okay, so this removed all the data from the tuning curves. I think there
%is a clear discrepancy in the size of the responses (changes over time and
%increases), so I think we should correct to the median value of the
%overall calculated peak values. Implementing this above. 


%hmmm....dont seem  to really see any kind of coherent difference. third
%mouse shows some changes, but these are very general and broad.
%Confusing...

for bigInd = 1:numFiles
    
end



