function [output] = functionWavePropertyExtraction(averageWaveData);

%hardcoded values. hopefully these shouldnt change.
trodesFS = 30000;
microSec = 1000000;

%want to get finer resolution (microseconds). Therefore, need to pull a
%linear interpolation.

%pull the samples to set up for interpolation
trodesSamples = (1:1:40)/trodesFS;
microSamples = 1/trodesFS:1/microSec:40/trodesFS;
interpWaves = zeros(size(microSamples,2),size(averageWaves,2));
%interpolates to 1 us resolution.
for waveCount = 1:size(averageWaves,2)
    interpWaves(:,waveCount) = interp1(trodesSamples,averageWaves(:,1),microSamples);
end
averageWaves = interpWaves;

%generates negative waves to pick up troughs
averageWavesNeg = averageWaves*-1;

%cell array for storage of data.
waveData = cell(size(averageWaves,2),3);

for waveCount = 1:size(averageWaves,2)
    %find peaks for signal in the correct orientation
    [pks locs] = findpeaks(averageWaves(:,waveCount));
    %stores things in an array: columns are different waves. row 1 is peak
    %size, row 2 is peak location, row 3 is peak half width (in
    %microseconds)
    peaksPos = zeros(3,size(pks,1));
    peaksPos(1,:) = pks';
    peaksPos(2,:) = locs';
    %for loop to find halfwidth. determines half height, calculates nearest
    %points around this half height.
    for peakCount = 1:size(pks,1) 
        halfPeakSize = peaksPos(1,peakCount)/2;
        leftMark = find(averageWaves(1:peaksPos(2,peakCount),waveCount)<halfPeakSize,1,'last');
        rightMark = find(averageWaves(peaksPos(2,peakCount):end,waveCount)<halfPeakSize,1,'first')+peaksPos(2,peakCount);
        peaksPos(3,peakCount) = rightMark-leftMark;
    end
    %find peaks for the negative image of the average wave. Does the same
    %operations.
    [pksneg locsneg] = findpeaks(averageWavesNeg(:,waveCount));
    %if statement is safeguard against waveforms with no negative peaks.
    if size(pksneg,1) > 0
        peaksNeg = zeros(3,size(pksneg,1));
        peaksNeg(1,:) = pksneg';
        peaksNeg(2,:) = locsneg';
        %for loop to find halfwidth. determines half height, calculates nearest
        %points around this half height.
        for peakCount = 1:size(pksneg,1) 
            halfPeakSize = peaksNeg(1,peakCount)/2;
            leftMark = find(averageWavesNeg(1:peaksNeg(2,peakCount),waveCount)<halfPeakSize,1,'last');
            rightMark = find(averageWavesNeg(peaksNeg(2,peakCount):end,waveCount)<halfPeakSize,1,'first')+peaksNeg(2,peakCount);
            peaksNeg(3,peakCount) = rightMark-leftMark;
        end
        waveData{waveCount,1} = peaksPos;
        waveData{waveCount,2} = peaksNeg;
        waveData{waveCount,3} = size(peaksPos,2)+size(peaksNeg,2); %pulls total number of peaks
    else
        %pulls total number of peaks
        waveData{waveCount,3} = size(peaksPos,2);
    end
    
end

%create structured array for storage
output = struct;
output.PeakInfo = waveData;

end