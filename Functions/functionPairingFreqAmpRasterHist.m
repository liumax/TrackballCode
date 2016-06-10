%This function is meant to generate average rasters and histograms for each
%frequency/amplitude pair. Also meant to calculate the overall response
%reliability to the cued event, dividing response reliability with the
%expected % of trials with a spike just by random chance. 

function [structureTarget] = functionPairingFreqAmpRasterHist(i, clusters,...
    structureTarget,histBinNum,histBinVector,histBin,truncatedNames,...
    master,soundData,masterStruct);

uniqueDBs = soundData.UniqueDBs;
uniqueFreqs = soundData.UniqueFreqs;
numDBs = size(uniqueDBs,1);
numFreqs = size(uniqueFreqs,1);

respReliab = cell(clusters,1);

for k = 1:clusters
    percentResponse = zeros(numFreqs,numDBs,1);
    for j=1:numFreqs;
        for l = 1:numDBs;
            aveFreqHolder = zeros(histBinNum,1);
            %finds all trials of specific frequency and amplitude
            trialNumHolder = intersect(find(master(:,2) == uniqueFreqs(j)),...
        find(master(:,3) == uniqueDBs(l)));
    %generates placeholder for responses per trial
            responseHolder = zeros(size(trialNumHolder,1),2);
            responseHolder(:,1) = trialNumHolder;
            %this goes through the individual trials and finds if there
            %are any responses during the tone period
            for m = 1:soundData.ToneReps;
                indivResp = size(find(structureTarget.Rasters{k}(:,2)>0 & ...
        structureTarget.Rasters{k}(:,2)<soundData.ToneDur &...
        structureTarget.Rasters{k}(:,1)==responseHolder(m)),1);
                responseHolder(m,2) = indivResp;
            end
            %uses these values to calculate response reliability,
            %stores this as a percentage
            percentResponse(j,l) = (soundData.ToneReps - size(find(responseHolder(:,2) == 0),1))...
                /soundData.ToneReps;
            %this pulls all times from all trials for the given dB and
            %frequency, and puts them all together 
            spikeTimeHolder = structureTarget.Rasters{k}...
                (structureTarget.Rasters{k}(:,3) == ...
            soundData.UniqueFreqs(j) & ...
            structureTarget.Rasters{k}(:,4) ...
            == soundData.UniqueDBs(l),2);
        %converts this data into the form of a histogram for storage
            [counts centers] = hist(spikeTimeHolder,histBinVector);
            countSize = size(counts);
            centerSize = size(centers);
            %this is to compensate for issues with the output being a
            %column/vector
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            freqDBHist{j,l} = [counts'*(1/histBin)*(1/soundData.ToneReps),centers'];%This should make it such taht units are in Hz/trial overall. That is, the average response per trial.
            aveFreqHolder(:,1) = aveFreqHolder(:,1) + freqDBHist{j,l}(:,1)/size(soundData.UniqueDBs,1); %160529 adjusted so that this doesnt produce overinflated values due to repetitions of multiple dbs.
            counts = [];
            centers = [];
        end
        averageFreqResp(:,j) = aveFreqHolder;
        aveFreqHolder = [];
    end
    masterFreqDBHist{k} = freqDBHist;
    expectedResponse = masterStruct.(truncatedNames{i}).AverageFiringRates(k)*soundData.ToneDur;
    if expectedResponse > 1
        expectedResponse = 1;
    end
    respReliab{k} = percentResponse/expectedResponse;
    aveFreqRespMaster{k} = averageFreqResp;
end

 %calculates responses binned per tone frequency
processRaster = cell(clusters,1);

%this pulls out correct frequency and correct timing, then finds which
%indices are correct for both.

masterToneRaster = structureTarget.Rasters;
for j = 1:clusters
    processRaster{j} = zeros(length(uniqueDBs),length(uniqueFreqs));
    for k = 1:length(uniqueFreqs)
        for m = 1:length(uniqueDBs)
        corrFreq = find(masterToneRaster{j}(:,3) == uniqueFreqs(k));
        corrDBs = find(masterToneRaster{j}(:,4) == uniqueDBs(m));
        corrTime = find(masterToneRaster{j}(:,2) > 0 & masterToneRaster{j}(:,2) < soundData.ToneDur);
        processRaster{j}(m,k) = length(intersect(intersect(corrFreq,corrTime),corrDBs));  
        end
    end
end

structureTarget.FreqDBSpecificHist = masterFreqDBHist;
structureTarget.ResponseReliability = respReliab;
structureTarget.AverageFrequencyHistogram = aveFreqRespMaster;
structureTarget.FrequencyResponse = processRaster;

end