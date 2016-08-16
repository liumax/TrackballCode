%This function is meant to generate average rasters and histograms for each
%frequency/amplitude pair. Also meant to calculate the overall response
%reliability to the cued event, dividing response reliability with the
%expected % of trials with a spike just by random chance. 

function [matclustStruct] = functionFreqAmpRasterHist(i, clusterSizer,...
    matclustStruct,histBinNum,histBinVector,histBin,truncatedNames);

respReliab = cell(clusterSizer,1);

for k = 1:clusterSizer
    percentResponse = zeros(size(matclustStruct.UniqueFreqs,1),size(matclustStruct.UniqueDBs,1));
    for j=1:size(matclustStruct.UniqueFreqs,1);
        for l = 1:size(matclustStruct.UniqueDBs,1);
            aveFreqHolder = zeros(histBinNum,1);
            %finds all trials of specific frequency and amplitude
            trialNumHolder = intersect(find(matclustStruct.Frequencies == matclustStruct.UniqueFreqs(j)),...
        find(matclustStruct.dBs == matclustStruct.UniqueDBs(l)));
    %generates placeholder for responses per trial
            responseHolder = zeros(size(trialNumHolder,1),2);
            responseHolder(:,1) = trialNumHolder;
            %this goes through the individual trials and finds if there
            %are any responses during the tone period
            for m = 1:matclustStruct.ToneReps;
                indivResp = size(find(matclustStruct.(truncatedNames{i}).Rasters{k}(:,2)>0 & ...
        matclustStruct.(truncatedNames{i}).Rasters{k}(:,2)<matclustStruct.ToneDur &...
        matclustStruct.(truncatedNames{i}).Rasters{k}(:,1)==responseHolder(m)),1);
                responseHolder(m,2) = indivResp;
            end
            %uses these values to calculate response reliability,
            %stores this as a percentage
            percentResponse(j,l) = (matclustStruct.ToneReps - size(find(responseHolder(:,2) == 0),1))...
                /matclustStruct.ToneReps;
            %this pulls all times from all trials for the given dB and
            %frequency, and puts them all together 
            spikeTimeHolder = matclustStruct.(truncatedNames{i}).Rasters{k}...
                (matclustStruct.(truncatedNames{i}).Rasters{k}(:,3) == ...
            matclustStruct.UniqueFreqs(j) & ...
            matclustStruct.(truncatedNames{i}).Rasters{k}(:,4) ...
            == matclustStruct.UniqueDBs(l),2);
        %converts this data into the form of a histogram for storage
            [counts centers] = hist(spikeTimeHolder,histBinVector);
            countSize = size(counts);
            centerSize = size(centers);
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            freqDBHist{j,l} = [counts'*(1/histBin)*(1/matclustStruct.ToneReps),centers'];%This should make it such taht units are in Hz/trial overall. That is, the average response per trial.
            aveFreqHolder(:,1) = aveFreqHolder(:,1) + freqDBHist{j,l}(:,1)/size(matclustStruct.UniqueDBs,1); %160529 adjusted so that this doesnt produce overinflated values due to repetitions of multiple dbs.
            counts = [];
            centers = [];
        end
        averageFreqResp(:,j) = aveFreqHolder;
        aveFreqHolder = [];
    end
    masterFreqDBHist{k} = freqDBHist;
    expectedResponse = matclustStruct.(truncatedNames{i}).AverageFiringRate(k)*matclustStruct.ToneDur;
%     respReliab{k} = percentResponse/expectedResponse;%removed this expected response division because it generates uninterpretable outputs.
    respReliab{k} = percentResponse;
    aveFreqRespMaster{k} = averageFreqResp;
end

 %calculates responses binned per tone frequency
processRaster = cell(clusterSizer,1);

%this pulls out correct frequency and correct timing, then finds which
%indices are correct for both.
uniqueDBs = matclustStruct.UniqueDBs;
uniqueFreqs = matclustStruct.UniqueFreqs;
masterToneRaster = matclustStruct.(truncatedNames{i}).Rasters;
for j = 1:clusterSizer
    processRaster{j} = zeros(length(uniqueDBs),length(uniqueFreqs));
    for k = 1:length(uniqueFreqs)
        for m = 1:length(uniqueDBs)
        corrFreq = find(masterToneRaster{j}(:,3) == uniqueFreqs(k));
        corrDBs = find(masterToneRaster{j}(:,4) == uniqueDBs(m));
        corrTime = find(masterToneRaster{j}(:,2) > 0 & masterToneRaster{j}(:,2) < matclustStruct.ToneDur);
        processRaster{j}(m,k) = length(intersect(intersect(corrFreq,corrTime),corrDBs));  
        end
    end
end

matclustStruct.(truncatedNames{i}).FreqDBSpecificHist = masterFreqDBHist;
matclustStruct.(truncatedNames{i}).ResponseReliability = respReliab;
matclustStruct.(truncatedNames{i}).AverageFrequencyHistogram = aveFreqRespMaster;
matclustStruct.(truncatedNames{i}).FrequencyResponse = processRaster;

end