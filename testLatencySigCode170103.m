


%testing with some really obvious cells first:
cd E:\170103TestTuningSignificanceAndTimingDetection\161216_ML161110C_R17_3200_fastTuni
fileName = '161216_ML161110C_R17_3200_fastTuning'
[s] = analysisBasicTuning(fileName);

%lets look at 7 cluster 1 first, since this cell has background firing, and
%a really nice response. 

tester = s.nt7cluster1.AllRasters(:,1);

% %lets look at 7 cluster 2 , no background firing
% 
% tester = s.nt7cluster2.AllRasters(:,1);

histBin = 0.001; %set this as a default

rasterWindow = s.Parameters.RasterWindow; %pull the appropriate raster window from params
periodBaseline = [rasterWindow(1) 0];
periodTone = [0 0.1];
periodPostTone = [0.1 0.2]; 
binsTone = (periodTone(2)-periodTone(1))/histBin;
binsPostTone = (periodPostTone(2)-periodPostTone(1))/histBin;

histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];
toneStartBin = find(histBinVector > 0,1,'first');

overallHist = hist(tester,histBinVector);

%plot to see
figure
plot(histBinVector,overallHist)
hold on
plot(histBinVector,overallHist,'r*')

%find percentile for 0.001 rank (0.1%)
percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);


%looks like there are about 4-5 bins that are above the baseline
pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
pHolderPostTone = find(overallHist(toneStartBin+binsTone:toneStartBin+binsTone+binsPostTone-1)>percentLimit);


%looks good with percentile!

%try with another dataset
tester = s.nt7cluster2.AllRasters(:,1);
%still looks good! how about a dataset that would be shit with such tiny
%bins

%2 cluster 1 is a bit more smeared out
tester = s.nt2cluster1.AllRasters(:,1);
%looks good too


%Moving forward with percentiles, since rank sum seems to break even with
%fairly low key normal distributions, which suggests that either I'm using
%it wrong, or i need some other test. 

%now look at latency.

%First, I need to determine whether there are multiple bins that are
%"significant"
respToneLength = length(pHolderTone);
%then I need to find if these values are consecutive
respToneDiff = diff(pHolderTone);
%find first of consecutive string
respToneFirst = find(respToneDiff == 1, 1, 'first');

if respToneLength > 1 && ~isempty(respToneFirst)
    percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
    latSwitch = 0; %switch to get out of while loop
    latInd = 1; %index that increments to move backwards in time
    while latSwitch == 0
        respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
        if respVal > percentLimitLow
            latInd = latInd + 1;
        elseif respVal < percentLimitLow
            latInd = latInd - 1;
            latSwitch = 1;
        end
        if latInd >= pHolderTone(respToneFirst)
            latSwitch = 1;
        end
    end
    if latInd <= respToneFirst
        toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
    else
        toneRespLat = 0;
    end
end

%now lets find peak response, which is fairly simple
toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset


%Now lets try for different frequencies/dBs

uniqueFreqs = s.SoundData.UniqueFrequencies;
uniqueDBs = s.SoundData.UniqueDBs;

numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

tester = s.nt2cluster1.AllRasters(:,1);
overallHist = hist(tester,histBinVector);
percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt2cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%Now try with a more sharply responding cell

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

% tester = s.nt7cluster1.AllRasters(:,1);
% overallHist = hist(tester,histBinVector);
% percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt7cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end


%for 7 unit 1, needed to smooth over 3 bins, otherwise responses were too
%limited to form contiguous responses...

%Try 7 unit 2

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

% tester = s.nt7cluster1.AllRasters(:,1);
% overallHist = hist(tester,histBinVector);
% percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt7cluster2.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end


%lets try with separate dataset:
clear
cd E:\170103TestTuningSignificanceAndTimingDetection\161216_ML161110C_R17_3000_fastTuni
fileName = '161216_ML161110C_R17_3000_fastTuning'
[s] = analysisBasicTuning(fileName);

tester = s.nt2cluster1.AllRasters(:,1);

histBin = 0.001; %set this as a default

rasterWindow = s.Parameters.RasterWindow; %pull the appropriate raster window from params
periodBaseline = [rasterWindow(1) 0];
periodTone = [0 0.1];
periodPostTone = [0.1 0.3]; 
binsTone = (periodTone(2)-periodTone(1))/histBin;
binsPostTone = (periodPostTone(2)-periodPostTone(1))/histBin;

histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];
toneStartBin = find(histBinVector > 0,1,'first');

overallHist = hist(tester,histBinVector);

uniqueFreqs = s.SoundData.UniqueFrequencies;
uniqueDBs = s.SoundData.UniqueDBs;

numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt2cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end
%works fairly well, maybe some of the responses are exaggerated with this
%method

%try 5 cluster 1 (should have no response, but possible noise)

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt5cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end
%works fairly well, picks up a bit of noise

%try 5 cluster 3

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt5cluster3.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%as expected, picks up a single noise point

%try 6 cluster 1 (really sharply aligned response)

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt6cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%again, works fairly well, does pick up some noise

%try 8 cluster 1, which is a not super well tuned neuron

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt8cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%again, works okay, picks up one noise

%try 7 cluster 2, should just be noise, but has some stuff that might get
%picked up

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt7cluster2.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%as expected, picked up the noise, but otherwise clean :D

%lets do one final dataset
clear
cd E:\170103TestTuningSignificanceAndTimingDetection\160909_ML160718B_R_2200_fastTuni
fileName = '160909_ML160718B_R_2200_fastTuning'
[s] = analysisBasicTuning(fileName);

tester = s.nt28cluster1.AllRasters(:,1);

histBin = 0.001; %set this as a default

rasterWindow = s.Parameters.RasterWindow; %pull the appropriate raster window from params
periodBaseline = [rasterWindow(1) 0];
periodTone = [0 0.1];
periodPostTone = [0.1 0.3]; 
binsTone = (periodTone(2)-periodTone(1))/histBin;
binsPostTone = (periodPostTone(2)-periodPostTone(1))/histBin;

histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];
toneStartBin = find(histBinVector > 0,1,'first');

overallHist = hist(tester,histBinVector);

uniqueFreqs = s.SoundData.UniqueFrequencies;
uniqueDBs = s.SoundData.UniqueDBs;

numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);


for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt28cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%kind of dirty, but reasonably picking things up

%try another unit, which should have no response

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);


for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt25cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end


%picks up multiple points....

%try another unit, 19 cluster 1, which has response

latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);


for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt19cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9)
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%oddly enough, misses a bunch of responses
%try again, but use percentile cutoff based from aggregate data


latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);

tester = s.nt19cluster1.AllRasters(:,1);
overallHist = hist(tester,histBinVector);
percentLimit = prctile(overallHist(1:toneStartBin-1),99.9)/(numFreqs*numDBs);

for i = 1:numFreqs
    for j = 1:numDBs
        tester = s.nt19cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
%         percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time
            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end
            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)));
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end
    end
end

%This method is a piece of shit. Do not use. 

%Stick with using specific preceding segments for baseline calculations.
%Will cause some losses, but hey, thats life. 

%With this, I think the latency, and specific peak response time and
%magnitude calculation system is more or less complete. 

%whats my next step? I think for tuning curve calculations, its just a
%matter of binning over the tone period, and trying to detect if there is a
%response outside of the tone period. 

%lets try this by extending the duration I look at
periodTone = [0 0.3];
periodPostTone = [0.1 0.2]; 
binsTone = (periodTone(2)-periodTone(1))/histBin;


latStore = zeros(numFreqs,numDBs);
peakStore = zeros(numFreqs,numDBs);
peakTimeStore = zeros(numFreqs,numDBs);


for i = 1:numFreqs
    for j = 1:numDBs

        tester = s.nt28cluster1.FreqDBRasters{i,j}(:,1);
        overallHist = hist(tester,histBinVector);
        overallHist = smooth(hist(tester,histBinVector),3);
        percentLimit = prctile(overallHist(1:toneStartBin-1),99.9);
        pHolderTone = find(overallHist(toneStartBin:toneStartBin+binsTone-1)>percentLimit); %-1 to make 100 bins
        respToneLength = length(pHolderTone);
        %then I need to find if these values are consecutive
        respToneDiff = diff(pHolderTone);
        %find first of consecutive string
        respToneFirst = find(respToneDiff == 1, 1, 'first');

        if respToneLength > 1 && ~isempty(respToneFirst)
            percentLimitLow = prctile(overallHist(1:toneStartBin-1),95); %set lower threshold bound for tracing back to start of response
            latSwitch = 0; %switch to get out of while loop
            latInd = 1; %index that increments to move backwards in time

            while latSwitch == 0
                respVal = overallHist(toneStartBin + pHolderTone(respToneFirst) - latInd); %calculates histogram value 
                if respVal > percentLimitLow
                    latInd = latInd + 1;
                elseif respVal <= percentLimitLow
                    latInd = latInd - 1;
                    latSwitch = 1;
                end
                if latInd >= pHolderTone(respToneFirst)
                    latSwitch = 1;
                end
            end

            if latInd <= pHolderTone(respToneFirst)
                toneRespLat = histBinVector(toneStartBin + pHolderTone(respToneFirst)-latInd);
            else
                toneRespLat = 0;
            end
            %now lets find peak response, which is fairly simple
            toneRespPeak = max(overallHist(toneStartBin + pHolderTone(1):toneStartBin + pHolderTone(end)-1)); %NEED TO HAVE -1 here!!!
            toneRespPeakTime = histBinVector(find(overallHist == toneRespPeak,1,'first')); %gives time of peak response in ms relative to tone onset
            latStore(i,j) = toneRespLat;
            peakStore(i,j) = toneRespPeak;
            peakTimeStore(i,j) = toneRespPeakTime;  
        end

    end
end

%This fix seems to more or less work. still dirty, getting some noise. Need
%to figure out how I can eliminate this kind of stuff
