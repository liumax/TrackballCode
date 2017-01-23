

%NOTE: spikeName and ttlName must be in quotations.
function [masterStruct] = functionPairingToneAnalysis(masterStruct,desigNames,...
    spikeName,soundName,params); 
%pulls out things into variables to make them easier to call.
rasterWindow = params.rasterWindow;
toneWindow = params.toneWindow;
genWindow = params.genWindow;
latBin = params.latBin;
percentCutoff = params.percentCutoff;
baselineCutoff = params.baseCutoff;
histBin = params.histBin;
baselineBin = params.baselineBin;
calcWindow = params.calcWindow;
firstSpikeWindow = params.firstSpikeWindow;

soundData = masterStruct.SoundData.(soundName);
%pulls out freq and db data for easy access.
targetSet = soundData.TargetSet;
controlSet = soundData.ControlSet;
%pulls frequencies.
targetFreq = targetSet(1);
controlFreq = controlSet(1);
%pulls master array.
master = masterStruct.SoundData.(soundName).MasterArray;
%finds all TTL times and feeds them into the master array.
TTLTimes = masterStruct.TTLs.(soundName);
master(:,1) = TTLTimes;
%divide TTLs into two arrays: one for control, one for target.
TTLDesig = cell(2,1);
TTLDesig{1} = TTLTimes(soundData.Frequencies == targetFreq);
TTLDesig{2} = TTLTimes(soundData.Frequencies == controlFreq);
TTLNames = {'Target','Control'};

%pull tone repetitions
toneReps = soundData.ToneReps;

%converts raster window from ratio to actual time in seconds.
rasterWindow = rasterWindow*soundData.ToneDur;
toneWindow = toneWindow * soundData.ToneDur;
genWindow = genWindow * soundData.ToneDur;
calcWindow = calcWindow*soundData.ToneDur;
baselineBin = baselineBin*soundData.ToneDur;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];

numUnits = size(s.DesignationName,2);

s.BaselineSpikes.(spikeName) = zeros(numUnits,1);
s.SignificantSpikes.(spikeName) = zeros(numUnits,1);

for i = 1:numUnits
    %prep individual sets of spikes:
    spikeTimes = s.(desigNames{i}).(spikeName);
    %split into target and control data
    for j = 1:length(TTLDesig)
        disp(strcat('Calculating for:',TTLNames{j}))
        fieldName = strcat(soundName,TTLNames{j},'Analysis');
        %rasterize spike data
        [rasters] = functionBasicRaster(spikeTimes,TTLDesig{j},rasterWindow);
        %pull individual histograms
        indivHist = zeros(length(histBinVector),toneReps);
        averageSpike = zeros(toneReps,1);
        for k = 1:toneReps
            indivHist(:,k) = hist(rasters(rasters(:,2) == k,1),histBinVector)/histBin;
            averageSpike(k) = size(find(rasters(:,2) == k & rasters(:,1) > baselineBin(1) & rasters(:,1) < baselineBin(2)),1);
        end
        averageRate = mean(averageSpike/(baselineBin(2)-baselineBin(1)));
        averageSTD = std(averageSpike/(baselineBin(2)-baselineBin(1)));
        %make compensatory binning relative to average rate.
        compBinTone = averageRate*(toneWindow(2)-toneWindow(1));
        compBinGen = averageRate*(genWindow(2)-genWindow(1));
        %perform histogram of all trials
        [histCounts histCenters] = hist(rasters(:,1),histBinVector);
        fullHistData = histCounts'/toneReps/histBin;
        %calculate significant response for combined histogram of all responses
        %to all sounds.
        [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,TTLDesig{j},toneReps);

        disp(strcat('Baseline Spikes:',num2str(generalResponseHist.SpikeNumber),' Unit:',(desigNames{i})))
        s.BaselineSpikes.(spikeName)(i) = generalResponseHist.SpikeNumber;
        if generalResponseHist.Warning == 0 & generalResponseHist.SigSpike == 1
            s.SignificantSpikes.(spikeName)(i) = 1;
        end
        [LatPeakBinOut] = functionLatPeakBinCalculation(toneWindow,genWindow,rasterWindow,rasters,...
                toneReps,2,[1:1:toneReps],latBin,histBin,percentCutoff,baselineCutoff);
            
        s.(desigNames{i}).(fieldName).Rasters = rasters;
        s.(desigNames{i}).(fieldName).Histograms = fullHistData;
        s.(desigNames{i}).(fieldName).IndivHistograms = indivHist;

        s.(desigNames{i}).(fieldName).AverageRate = averageRate;
        s.(desigNames{i}).(fieldName).AverageSTD = averageSTD;

        s.(desigNames{i}).(fieldName).AllHistogramSig = generalResponseHist;
        s.(desigNames{i}).(fieldName).CompBinTone = compBinTone;
        s.(desigNames{i}).(fieldName).CompBinGen = compBinGen;

        s.(desigNames{i}).(fieldName).LatBinPeakCalcs = LatPeakBinOut;

        s.(desigNames{i}).(fieldName).BinToneComp = LatPeakBinOut.BinnedSpikesTone - compBinTone;
        s.(desigNames{i}).(fieldName).BinGenComp = LatPeakBinOut.BinnedSpikesGen - compBinGen;
        s.(desigNames{i}).(fieldName).PeakGenComp = LatPeakBinOut.PeakRespGen - averageRate;


        s.(desigNames{i}).(fieldName).HistBinVector = histBinVector;
    end
end

end