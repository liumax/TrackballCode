function [overOut,indivOut,s,masterData,masterHeader,masterInd] = functionTuningDataExtraction(numUnits,numDBs,numFreqs,uniqueFreqs,s,masterData,masterHeader,masterInd,histBinVector,totalTrialNum,master,desigNames,calcWindow,histBinNum,whiteStatus);
overOut.PosWidths = zeros(numDBs,numUnits,3); %store width of response for positive responses here. 
overOut.NegWidths = zeros(numDBs,numUnits,3); %store width of response for negative responses here. 
overOut.PosTots = zeros(numDBs,numUnits,3); %store total number of significant responses here
overOut.NegTots = zeros(numDBs,numUnits,3);  %store total number of significant responses here
overOut.PosCont = zeros(numDBs,numUnits,3); %store whether response is contiguous or not here
overOut.NegCont = zeros(numDBs,numUnits,3); %store whether response is contiguous or not here
overOut.PosEdgeWarn = zeros(numDBs,numUnits,3); %store whether significant responses abut an edge
overOut.NegEdgeWarn = zeros(numDBs,numUnits,3); %store whether significant responses abut an edge
overOut.PosGaussWidth = zeros(numDBs,numUnits,3); % will perform gaussian fit on binned spikes if there is sufficient significant values. Store half-peak width
overOut.RespWidthPos = zeros(numFreqs,numDBs,numUnits,2); %width of firing rate, positive
widthStore = zeros(numFreqs,numDBs,numUnits);
bigWidth = [];
for i = 1:numUnits
    masterHolder = masterInd;
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    alignTimes = master(:,1);
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Distance(1,1)/s.Parameters.trodesFS:s.Parameters.SpeedFiringBins:s.RotaryData.Distance(end,1)/s.Parameters.trodesFS]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    
    %calculates rasters based on spike information. 
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.RasterWindow);
    rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
    fullRasterData = rasters;
    
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    averageSpikeHolder = zeros(totalTrialNum,1);
    for k = 1:totalTrialNum
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > s.Parameters.BaseWindow(1) & rasters(:,1) < 0),1);
    end
    averageRate = mean(averageSpikeHolder/(-s.Parameters.BaseWindow(1)));
    averageSTD = std(averageSpikeHolder/(-s.Parameters.BaseWindow(1)));
    averageSTE = averageSTD/(sqrt(totalTrialNum-1));
    
    %store average rate into master
    masterData(i,masterHolder) = averageRate;
    masterHeader{masterHolder} = 'Average Rate';
    masterHolder = masterHolder + 1;
    
    %now pull waveform data and isi data to categorize cells
    try
        waveForms = s.(desigNames{i}).AverageWaveForms;
        
        %find size of waveforms
        waveSize = size(waveForms);

        %find maximum waveform by max peak size
        maxWave = max(waveForms);
        [maxVal maxInd] = max(maxWave);

        %chose the big wave, interpolate to fine degree
        chosenWave = waveForms(:,maxInd);
        interpVect = [1:0.1:40];
        interpWave = interp1(1:40,chosenWave,interpVect,'spline');
    catch
        waveForms = s.(desigNames{i}).medianWave*-1;
        
        %find size of waveforms
        waveSize = size(waveForms);

        %find maximum waveform by max peak size
        maxWave = max(waveForms');
        [maxVal maxInd] = max(maxWave);

        %chose the big wave, interpolate to fine degree
        chosenWave = waveForms(maxInd,:);
        interpVect = [1:0.1:length(chosenWave)];
        interpWave = interp1(1:length(chosenWave),chosenWave,interpVect,'spline');
    end
    
    
    %now we need to find the peak. Find this starting at point 10. 

    [pkVal pkInd] = max(interpWave(100:end));
    pkInd = pkInd + 100 - 1;
    %now we need to find the minimum following the peak

    [troughVal troughInd] = min(interpWave(pkInd:end));
    troughInd = troughInd + pkInd - 1;

    peakTrough = (troughInd - pkInd)/300000;
    
    %find ISIs
    isiTimes = diff(spikeTimes);
    isiCov = std(isiTimes)/mean(isiTimes);
    
    masterData(i,masterHolder) = peakTrough;
    masterHeader{masterHolder} = 'PeakTrough(ms)';
    masterHolder = masterHolder + 1;
    
    masterData(i,masterHolder) = isiCov;
    masterHeader{masterHolder} = 'isiCov';
    masterHolder = masterHolder + 1;
    
    
    if peakTrough < s.Parameters.PVLim(1) & isiCov > s.Parameters.ChatLim
        masterData(i,masterHolder) = 1; %pv cell
    elseif peakTrough > s.Parameters.PVLim(2) & isiCov < s.Parameters.ChatLim & averageRate > 2
        masterData(i,masterHolder) = 2; %ChAT Cell
    elseif peakTrough > s.Parameters.PVLim(2) & isiCov > s.Parameters.ChatLim
        masterData(i,masterHolder) = 0; %MSN
    else
        masterData(i,masterHolder) = NaN; %label as unknown
    end
    masterHeader{masterHolder} = 'CellType';
    masterHolder = masterHolder + 1;
    
    
    %calculates histograms of every single trial. 
    fullHistHolder = zeros(length(histBinVector),totalTrialNum);
    for k =1:totalTrialNum
        if size(rasters(rasters(:,2) == k,:),1) > 0
            [counts centers] = hist(rasters(rasters(:,2) == k,1),histBinVector);
            fullHistHolder(:,k) = counts;
        end
    end
    
    %calculate standard deviation for each bin
    histSTE = (std(fullHistHolder'))';
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/totalTrialNum/s.Parameters.histBin;
    
    %calculate significant response for combined histogram of all responses
    %to all sounds.    
%     [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,length(master));
    [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,totalTrialNum,...
        s.Parameters.minSpikes,s.Parameters.latBin,[s.Parameters.BaseWindow(1),0],s.Parameters.zLimit,s.Parameters.minSigSpikes,s.Parameters.SigSmoothWindow);
    
    disp(strcat('Baseline Spikes:',num2str(generalResponseHist.SpikeNumber),' Unit:',(desigNames{i})))
    s.BaselineSpikes = generalResponseHist.SpikeNumber;
    if generalResponseHist.Warning == 0 & generalResponseHist.SigSpikePos == 1 
        s.SignificantSpikes(i) = 1;
    end
    
    masterData(i,masterHolder) = generalResponseHist.MeanBaseline;
    masterHeader{masterHolder} = 'BaselineRate';
    masterHolder = masterHolder + 1;
    
    masterData(i,masterHolder) = generalResponseHist.SigSpikePos;
    masterHeader{masterHolder} = 'PosSigGenHist';
    masterHolder = masterHolder + 1;
    
    masterData(i,masterHolder) = generalResponseHist.SigSpikeNeg;
    masterHeader{masterHolder} = 'NegSigGenHist';
    masterHolder = masterHolder + 1;
    
    %allocates empty array.
    organizedHist = zeros(numFreqs,numDBs,length(histBinVector));
    organizedRasters = cell(numFreqs,numDBs);
    responseHistHolder = cell(numFreqs,numDBs);
    histErr = zeros(numFreqs,numDBs,length(histBinVector));
    latPeakBinStore = cell(numFreqs,numDBs);
    latStore = zeros(numFreqs,numDBs);
    peakStoreGen = zeros(numFreqs,numDBs);
    binStoreTone = zeros(numFreqs,numDBs);
    binStoreGen = zeros(numFreqs,numDBs);
    probStoreTone = zeros(numFreqs,numDBs);
    probStoreGen = zeros(numFreqs,numDBs);
    freqSpecHist = zeros(numFreqs,histBinNum,1);
    
    
    for k = 1:numFreqs
        subUniqueDB = unique(master(master(:,2) == uniqueFreqs(k),3));
        for l = 1:numDBs
            targetTrials = master(master(:,2) == uniqueFreqs(k) & master(:,3) == subUniqueDB(l),4); %finds the trial number for all trials of given frequency and amplitude
            findMatches = find(ismember(fullRasterData(:,2),targetTrials)); %uses these trial numbers to pull raster index
            targetRasters = fullRasterData(findMatches,:); %extracts target rasters from all rasters using the above index.
            organizedRasters{k,l} = targetRasters; %saves to organized rasters
            [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector); %calculates histogram with defined bin size
            organizedHist(k,l,:) = histCounts/length(targetTrials)/s.Parameters.histBin; %saves histogram
            specHist = fullHistHolder(:,targetTrials);
            histErr(k,l,:) = std(specHist,0,2)/sqrt(length(targetTrials));
            [latPeakBinOut] = functionLatPeakBinCalculation(s.Parameters.ToneWindow,s.Parameters.GenWindow,s.Parameters.RasterWindow,...
                targetRasters,length(targetTrials),2,targetTrials,s.Parameters.latBin,s.Parameters.PercentCutoff,s.Parameters.BaselineCutoff);
            latePeakBinStore{k,l} = latPeakBinOut;
            latStore(k,l) = latPeakBinOut.ResponseLatency;
            peakStoreFast(k,l) = latPeakBinOut.PeakRespFast;
            peakStoreTone(k,l) = latPeakBinOut.PeakRespTone;
            peakStoreGen(k,l) = latPeakBinOut.PeakRespGen;
            binStoreFast(k,l) = mean(latPeakBinOut.BinnedSpikesFast);
            binStoreTone(k,l) = mean(latPeakBinOut.BinnedSpikesTone);
            binStoreGen(k,l) = mean(latPeakBinOut.BinnedSpikesGen);
            probStoreTone(k,l) = latPeakBinOut.ProbSpikeTone;
            probStoreGen(k,l) = latPeakBinOut.ProbSpikeGen;
            binSigVals(k,l,:) = latPeakBinOut.BinSigVals;
            binDiff(k,l,:) = latPeakBinOut.BinDiff;
            latPeakBinOut = [];
%             [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,toneReps);
            [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes(targetTrials),length(targetTrials),...
        s.Parameters.minSpikes,s.Parameters.latBin,[s.Parameters.BaseWindow(1),0],s.Parameters.zLimit,s.Parameters.minSigSpikes,s.Parameters.SigSmoothWindow);
            responseHistHolder{k,l} = responseHist;
            if responseHist.SigSpikePos == 1%if have a significant response, check for width
                [widthOut] = functionResponseWidth(responseHist);
                if length(widthOut.Widths) > 0
                    widthStore(k,l,i) = widthOut.Widths(1);
                    widthLat(k,l,i) = widthOut.StartsEnds(1);
                    bigWidth{k,l,i} = widthOut;
                end
            end
        end
        if numDBs == 1
            freqSpecHist(k,:) = (squeeze(organizedHist(k,:,:)));
        elseif numDBs > 1
            freqSpecHist(k,:) = mean(squeeze(organizedHist(k,:,:)));
        end
        
    end
        disp(strcat('Raster and Histogram Extraction Complete: Unit ',num2str(i),' Out of ',num2str(numUnits)))
    indivOut.(desigNames{i}).AllRasters = fullRasterData;
    indivOut.(desigNames{i}).AllHistograms = fullHistData;
    indivOut.(desigNames{i}).IndividualHistograms = fullHistHolder; 
    indivOut.(desigNames{i}).HistogramStandardDeviation = histSTE;
    indivOut.(desigNames{i}).FreqDBRasters = organizedRasters;
    indivOut.(desigNames{i}).FreqDBHistograms = organizedHist;
    indivOut.(desigNames{i}).FreqDBHistogramErrors = histErr;
    indivOut.(desigNames{i}).FrequencyHistograms = freqSpecHist;
    indivOut.(desigNames{i}).AverageRate = averageRate;
    indivOut.(desigNames{i}).SessionFiring = sessionFiring;
    indivOut.(desigNames{i}).AverageSTD = averageSTD;
    indivOut.(desigNames{i}).AverageSTE = averageSTE;
    indivOut.(desigNames{i}).HistBinVector = histBinVector;
    indivOut.(desigNames{i}).AllHistogramSig = generalResponseHist;
    indivOut.(desigNames{i}).SpecHistogramSig = responseHistHolder;
    indivOut.(desigNames{i}).LatPeakBin = latePeakBinStore;
    indivOut.(desigNames{i}).LatencyMap = latStore;
    indivOut.(desigNames{i}).PeakMapGen = peakStoreGen;
    indivOut.(desigNames{i}).PeakMapTone = peakStoreTone;
    indivOut.(desigNames{i}).BinFast = binStoreFast;
    indivOut.(desigNames{i}).BinTone = binStoreTone;
    indivOut.(desigNames{i}).BinGen = binStoreGen;
    indivOut.(desigNames{i}).BinDiff = binDiff;
    indivOut.(desigNames{i}).ProbTone = probStoreTone;
    indivOut.(desigNames{i}).ProbGen = probStoreGen;
    indivOut.(desigNames{i}).BinSigVals = binSigVals;
%     indivOut.WidthData = widthOut;
    
    %store some information about width
    CFTrig = 0;
    BFTrig = 0;
    for j = 1:numDBs
        sigThresh = 0.01;
        %pull out sign of change, and significance value
        if whiteStatus == 1
            targetSig = squeeze(binSigVals(2:end,j,:));
            targetSign = squeeze(binDiff(2:end,j,:));
            targetBins = squeeze(binStoreFast(2:end,j,:));
            
        else
            targetSig = squeeze(binSigVals(:,j,:));
            targetSign = squeeze(binDiff(:,j,:));
            targetBins = squeeze(binStoreFast(:,j,:));
        end

        for kk = 1:size(targetSig,2)
            %first, deal with positives
            findPos =  find(targetSig(:,kk) < sigThresh & targetSign(:,kk) > 0);
            if findPos%in the event of positive and significant events
                
                if CFTrig == 0
                    sigValLength = length(find(targetSig(:,kk) < sigThresh));
                    if whiteStatus == 1
                        if sigValLength == 1
                            masterData(i,masterHolder) = uniqueFreqs(find(targetSig(:,kk)<sigThresh)+1);
                        else
                            masterData(i,masterHolder) = mean(uniqueFreqs(find(targetSig(:,kk)<sigThresh)+1));
                        end
                    else
                        if sigValLength == 1
                            masterData(i,masterHolder) = uniqueFreqs(find(targetSig(:,kk)<sigThresh));
                        else
                            masterData(i,masterHolder) = mean(uniqueFreqs(find(targetSig(:,kk)<sigThresh)));
                        end
                    end
                    masterHeader{masterHolder} = 'CF';
                    masterHolder = masterHolder + 1;
                    disp('CF Detected')
                    CFTrig = 1;
                end
                %determine number of significant values
                if find(findPos == 1) | find(findPos == length(targetSig)) %determine if anything on the edge of the tuning range
                    overOut.PosEdgeWarn(j,i,kk) = 1;
                end
                if length(findPos) > 1
                    try
                        x = [1:length(targetSig)];
                        y = targetBins;
                        f = fit(x',y,'gauss1'); %perform single gaussian fit
                        newVect = [1:0.1:length(targetSig)];
                        fVals = f(newVect);
                        %find peak
                        [peakVal peakInd] = max(fVals);

                        %find half peak width
                        firstBound = find(fVals(1:peakInd) - peakVal/2 > 0,1,'first');
                        lastBound = find(fVals(peakInd:end) - peakVal/2 > 0,1,'last');
                        overOut.PosGaussWidth(j,i,kk) = (lastBound + peakInd - firstBound)/10;
                    catch
                        overOut.PosGaussWidth(j,i,kk) = 0;
                    end
                end
                overOut.PosTots(j,i,kk) = length(findPos);
                diffFind = diff(findPos);
                if diffFind == 1
                    disp('All Consecutive!')
                    overOut.PosCont(j,i,kk) = 1;
                    %since all consecutive, width equals length of findPos
                    overOut.PosWidths(j,i,kk) = length(findPos);
                else
                    disp('Not Consecutive, Finding Widest')
                    overOut.PosCont(j,i,kk) = 0;
%                     findSpaces = length(find(diffFind ~= 1)); %determine how many different segments there are
                    findGaps = find(diffFind ~= 1);
                    findGaps(2:end+1) = findGaps;
                    findGaps(1) = 0;
                    findGaps(end+1) = length(findPos);
                    diffLengths = diff(findGaps);
                    disp('Widest Point Found')
                    overOut.PosWidths(j,i,kk) = max(diffLengths);
                end
                %find maximum value. This is BF. 
                if j == numDBs & BFTrig == 0
                    [M bfFind] = max(findPos);
                    if whiteStatus == 1
                        masterData(i,masterHolder) = uniqueFreqs(findPos(bfFind)+1);
                    else
                        masterData(i,masterHolder) = uniqueFreqs(findPos(bfFind));
                    end
                    masterHeader{masterHolder} = 'BF';
                    masterHolder = masterHolder + 1;
                    disp('BF Detected')
                    BFTrig = 1;
                end
            else
                disp('No Positives Found')
            end
            
            %now do negative
            sigThresh = 0.01;
            findNeg =  find(targetSig(:,kk) < sigThresh & targetSign(:,kk) < 0);
            if findNeg%in the event of positive and significant events
                if find(findNeg == 1) | find(findNeg == length(targetSig))
                    overOut.NegEdgeWarn(j,i,kk) = 1;
                end
                overOut.NegTots(j,i,kk) = length(findNeg);
                diffFind = diff(findNeg);
                if diffFind == 1
                    disp('All Consecutive!')
                    overOut.NegCont(j,i,kk) = 1;
                    %since all consecutive, width equals length of findNeg
                    overOut.NegWidths(j,i,kk) = length(findNeg);
                else
                    disp('Not Consecutive, Finding Widest')
                    overOut.NegCont(j,i,kk) = 0;
%                     findSpaces = length(find(diffFind ~= 1)); %determine how many different segments there are
                    findGaps = find(diffFind ~= 1);
                    findGaps(2:end+1) = findGaps;
                    findGaps(1) = 0;
                    findGaps(end+1) = length(findNeg);
                    diffLengths = diff(findGaps);
                    disp('Widest Point Found')
                    overOut.NegWidths(j,i,kk) = max(diffLengths);
                end
            else
                disp('No Negatives Found')
            end
        end
        if j == numDBs & CFTrig == 0
            masterData(i,masterHolder) = NaN;
            masterHeader{masterHolder} = 'CF';
            masterHolder = masterHolder + 1;
            disp('NO CF')
        end
        
        if j == numDBs & BFTrig == 0
            masterData(i,masterHolder) = NaN;
            masterHeader{masterHolder} = 'BF';
            masterHolder = masterHolder + 1;
            disp('NO BF')
        end
    end
    
    
end
overOut.WidthData = widthStore;
overOut.WidthLatData = widthLat;
% overOut.FullWidth = bigWidth;

masterInd = masterHolder;

end
