function [output] = functionWavePropertyExtraction(averageWaves);

%hardcoded values. hopefully these shouldnt change.
trodesFS = 30000;
microSec = 1000000;
minPeakWidth = 50; %minimum half-width for any peak, in microseconds.

%want to get finer resolution (microseconds). Therefore, need to pull a
%linear interpolation.

%pull the samples to set up for interpolation
trodesSamples = (1:1:40)/trodesFS;
microSamples = 1/trodesFS:1/microSec:40/trodesFS;
interpWaves = zeros(size(microSamples,2),size(averageWaves,2));
%interpolates to 1 us resolution.
for waveCount = 1:size(averageWaves,2)
    interpWaves(:,waveCount) = interp1(trodesSamples,averageWaves(:,waveCount),microSamples);
end
averageWaves = interpWaves;

%generates negative waves to pick up troughs
averageWavesNeg = averageWaves*-1;

%generate array to store wave information. first row is half-width, second
%row is peak-trough time
waveInfo = zeros(2,size(averageWaves,2));

for waveCount = 1:size(averageWaves,2)
    %find peaks for signal in the correct orientation
    [pks locs] = findpeaks(averageWaves(:,waveCount));
    %stores things in an array: columns are different waves. row 1 is peak
    %size, row 2 is peak location, row 3 is peak half width (in
    %microseconds)
    peaksPos = zeros(size(pks,1),3);
    peaksPos(:,1) = pks';
    peaksPos(:,2) = locs';
    %for loop to find halfwidth. determines half height, calculates nearest
    %points around this half height.
    for peakCount = 1:size(pks,1) 
        halfPeakSize = peaksPos(peakCount,1)/2;
        leftMark = find(averageWaves(1:peaksPos(peakCount,2),waveCount)<halfPeakSize,1,'last');
        rightMark = find(averageWaves(peaksPos(peakCount,2):end,waveCount)<halfPeakSize,1,'first')+peaksPos(peakCount,2);
        peaksPos(peakCount,3) = rightMark-leftMark;
    end
    
    %now I insert a quality control check. This will remove all peaks that
    %are too narrow (below variable # of microseconds) and removes them! 
    findNarrow = find(peaksPos(:,3) < minPeakWidth);
    peaksPos(findNarrow,:) = [];
    
    %find peaks for the negative image of the average wave. Does the same
    %operations.
    [pksneg locsneg] = findpeaks(averageWavesNeg(:,waveCount));
    %if statement is safeguard against waveforms with no negative peaks.
    if size(pksneg,1) > 0
        peaksNeg = zeros(size(pksneg,1),3);
        peaksNeg(:,1) = (pksneg');
        peaksNeg(:,2) = locsneg';
        %for loop to find halfwidth. determines half height, calculates nearest
        %points around this half height.
        for peakCount = 1:size(pksneg,1) 
            halfPeakSize = peaksNeg(peakCount,1)/2;
            leftMark = find(averageWavesNeg(1:peaksNeg(peakCount,2),waveCount)<halfPeakSize,1,'last');
            rightMark = find(averageWavesNeg(peaksNeg(peakCount,2):end,waveCount)<halfPeakSize,1,'first')+peaksNeg(peakCount,2);
            peaksNeg(peakCount,3) = rightMark-leftMark;
        end
        %quality control to remove super minor peaks.
        findNarrow = find(peaksNeg(:,3) < minPeakWidth);
        peaksNeg(findNarrow,:) = [];
        peaksNeg(:,1) = -peaksNeg(:,1);
        
        allPeaks = [peaksPos;peaksNeg];
        allPeaks = sortrows(allPeaks,2);

    else
        allPeaks = peaksPos;
    end
    
    peakNums = size(allPeaks,1);
    
    if peakNums == 1;
        waveInfo(:,waveCount) = 0;
    elseif peakNums == 2;
        waveInfo(1,waveCount) = allPeaks(1,3);
        waveInfo(2,waveCount) = allPeaks(2,2) - allPeaks(1,2);
    elseif peakNums == 3;
        waveInfo(1,waveCount) = allPeaks(2,3);
        waveInfo(2,waveCount) = allPeaks(3,2) - allPeaks(2,2);
    elseif peakNums > 3;
        figure
        plot(averageWaves(:,waveCount))
        title('More than three peaks detected');
    end
end

%create structured array for storage
output = waveInfo;

end