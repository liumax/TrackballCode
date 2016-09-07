


baselineBinSize = 0.001; %hist bin in seconds
% preToneDuration = [-0.1 0];
% preToneBinCenters = [preToneDuration(1)+histBin/2:histBin:preToneDuration(2)-histBin/2];

toneHistBins = [0.01 0.001];
baselineBins = [];
toneRespWindow = [0 0.200];
preToneDuration = [-0.2 0];
preToneBinCenters = [preToneDuration(1)+baselineBinSize/2:baselineBinSize:preToneDuration(2)-baselineBinSize/2];
spikeProbWindow = [0 0.05 0.1 0.15]; %window over which spike responses will count towards probability of response

zLimit = 3; %z score for response to count as significant
% postToneBinCenters() = [postToneDuration(1)+toneHistBins/2:toneHistBins:postToneDuration(2)-toneHistBins/2]

% Looking at this: E:\160805ML160721\fastTuning\160802_ML160721A_R_3300_IDfastTuni
%ONE DATA SET
load('160802_ML160721A_R_3300_IDfastTuningFullTuningAnalysis.mat') %loads analysis
load('160802_ML160721A_R_3300_IDfastTuning.mat') %loads sound file
%pulls raster data out
%fairly good tuned cell. good firing rate increase, but not very time
%locked
rasterData = matclustStruct.matclust_param_nt7.Rasters{1,2};
%responsive cell, but extremely sparse. Very few spikes. Individual spikes
%pull z score to 4.0

% rasterData = matclustStruct.matclust_param_nt6.Rasters{1};
%responsive cell with low baseline firing. Extremely time locked response. 
% rasterData = matclustStruct.matclust_param_nt10.Rasters{1};

% rasterData = matclustStruct.matclust_param_nt10.Rasters{1};
%very responsive cell, almost no latency, basically single spike.


% %SECOND DATA SET
% load('160802_ML160721A_R_3300_IDfastTuningFullTuningAnalysis.mat')
% load('160802_ML160721A_R_3300_IDfastTuning.mat') %loads sound file




%plots to check for appropriate response
figure
plot(rasterData(:,2),rasterData(:,1),'b.')
%pulls frequency data
uniqueFreqs = unique(soundData.Frequencies);
uniqueDBs = unique(soundData.dBs);
toneReps = soundData.ToneRepetitions;
totalTrialNum = toneReps*size(uniqueFreqs,1)*size(uniqueDBs,1);

%one alternative to finding basline firing rate: average across all trials,
%then find mean and std from individual 5ms bins.
baselineHolder = rasterData(rasterData(:,2)<0,:);
%make histogram of baseline holder to generate baseline from which we can
%calculate response properties. 
[counts centers] = hist(baselineHolder(:,2),preToneBinCenters);
baselineCounts = counts/baselineBinSize/totalTrialNum;
baselineMean = mean(baselineCounts);
baselineSTD = std(baselineCounts);

%this version of the code makes use of smooth to generate smoothed data
%curves of the appropriate bin sizes.
rasterHolder = cell(size(uniqueFreqs,1),size(uniqueDBs,1));
histHolder = cell(size(uniqueFreqs,1),size(uniqueDBs,1));
rawHistHolder = cell(size(uniqueFreqs,1),size(uniqueDBs,1)); %will report purely the number of spikes per bin
sigHolder = zeros(size(uniqueFreqs,1),size(uniqueDBs,1),size(toneHistBins,2),5); 
%sigHolder stores information about spiking responses. in the fourth
%dimension, 1 stores first spike, 2 stores number of spikes exceeding z
%limit, 3 stores peak spike value in z score, 4 is peak in spikes/sec
%calculates what the firing rate is for a single spike. fifth is to report
%that there is no continuous peak
singleSpikeFreq = 1/0.001/toneReps;

%% goes through data and extracts features.
for i = 1:size(uniqueFreqs,1)
    for j = 1:size(uniqueDBs,1)
        rasterHolder{i,j} = rasterData(rasterData(:,3) == uniqueFreqs(i) & rasterData(:,4) == uniqueDBs(j) & rasterData(:,2) > toneRespWindow(1) & rasterData(:,2) <= toneRespWindow(2),:);
        histHolder{i,j} = cell(length(toneHistBins),1);
        rawHistHolder{i,j} = cell(length(toneHistBins),1);
        %generate histogram using 1ms bins. this is hardcoded.
        toneBinCenters = [toneRespWindow(1)+0.001/2:0.001:toneRespWindow(2)-0.001/2];
        [counts centers] = hist(rasterHolder{i,j}(:,2),toneBinCenters);
        counts = counts/0.001/toneReps; %this normalizes so it is calculated per trial
        for k = 1:size(toneHistBins,2)
            smoothSpan = toneHistBins(k)/0.001;
            if smoothSpan == 1;
                smoothedToneCounts = counts;
            else
                smoothedToneCounts = smooth(counts,smoothSpan);
            end
            rawHistHolder{i,j}{k} = smoothedToneCounts;

            histHolder{i,j}{k} = (smoothedToneCounts-baselineMean)/baselineSTD;%zscores the firing rate
            singleSpikeZ = (singleSpikeFreq-baselineMean)/baselineSTD;
            testHistHolder{i,j}{k} = zscore(smoothedToneCounts);
            if ~isempty(find(histHolder{i,j}{k} > zLimit))
                firstSig = find(histHolder{i,j}{k} > zLimit,1,'first');
%                 sigHolder(i,j,k,2) = length(find(histHolder{i,j}{k} > zLimit));
                findSigSpikes = find(histHolder{i,j}{k} > zLimit);
                sigSpikesDiff = diff(findSigSpikes);%finds continuous stretches of significant values
                if length(findSigSpikes) == 1
                    peakSize = histHolder{i,j}{k}(firstSig);
                    if peakSize == singleSpikeZ
                        sigHolder(i,j,k,1) = 0;
                        sigHolder(i,j,k,2) = 0;
                        sigHolder(i,j,k,3) = 0;
                        sigHolder(i,j,k,4) = 0;
                    else
                        sigHolder(i,j,k,2) = 1;
                        sigHolder(i,j,k,3) = peakSize;
                        sigHolder(i,j,k,4) = rawHistHolder{i,j}{k}(firstSig);
                    end
                elseif length(findSigSpikes) > 1 & ~isempty(find(sigSpikesDiff==1))
                    
                    firstSpikeTime = find(sigSpikesDiff==1,1,'first');%finds the first 1, which should be the beginning of the first true peak
                    firstSpikeTime = findSigSpikes(firstSpikeTime);
                    peakSize = max(histHolder{i,j}{k});
                    sigHolder(i,j,k,1) = firstSpikeTime;
                    sigHolder(i,j,k,2) = length(findSigSpikes);
                    sigHolder(i,j,k,3) = peakSize;
                    sigHolder(i,j,k,4) = max(rawHistHolder{i,j}{k});
                elseif length(findSigSpikes) > 1 & isempty(find(sigSpikesDiff==1))
                    sigHolder(i,j,k,1) = firstSig;
                    sigHolder(i,j,k,2) = length(findSigSpikes);
                    sigHolder(i,j,k,3) = peakSize;
                    sigHolder(i,j,k,4) = max(rawHistHolder{i,j}{k});
                    sigHolder(i,j,k,5) = 1;
                end

            end
        end
    end
end

%% some things to visualize results
% sigHolder(:,:,1,1)
% sigHolder(:,:,1,2)
% 
% sigHolder(:,:,2,1)
% sigHolder(:,:,2,2)
% 
% sigHolder(:,:,3,1)
% sigHolder(:,:,3,2)
% 
% sigHolder(:,:,4,1)
% sigHolder(:,:,4,2)
% 
% figure
% hold on
% plot(histHolder{7,1}{1,1})
% plot(histHolder{7,1}{2,1},'r')
% plot(histHolder{7,1}{3,1},'g')
% plot(histHolder{7,1}{4,1},'c')
% plot(histHolder{7,1}{4,1},'r*')
% figure
% hold on
% plot(rawHistHolder{7,1}{1,1})
% plot(rawHistHolder{7,1}{2,1},'r')
% plot(rawHistHolder{7,1}{3,1},'g')
% plot(rawHistHolder{7,1}{4,1},'c')
% plot(rawHistHolder{7,1}{4,1},'r*')

%% looking at average first spike timing, also probability of first spike

firstSpikeTimes = zeros(totalTrialNum,3);

for i = 1:totalTrialNum
    if ~isempty(find(rasterData(:,1) == i & rasterData(:,2) >0,1,'first'))
        findFirstSpike = find(rasterData(:,1) == i & rasterData(:,2) >0,1,'first');
        firstSpikeTimes(i,1) = rasterData(findFirstSpike,2);
        firstSpikeTimes(i,2) = rasterData(findFirstSpike,3);
        firstSpikeTimes(i,3) = rasterData(findFirstSpike,4);
    end
end

firstSpikeTimes(firstSpikeTimes(:,3) == 0,:) = [];

firstSpikeStats = zeros(length(uniqueFreqs),length(uniqueDBs),4);
%finds timing of the first spike in three time bins
for i = 1:size(uniqueFreqs,1)
    for j = 1:size(uniqueDBs,1)
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(i) & firstSpikeTimes(:,3) == uniqueDBs(j),1))
            firstSpikeStats(i,j,1) = mean(firstSpikeTimes(find(firstSpikeTimes(:,2) == uniqueFreqs(i) & firstSpikeTimes(:,3) == uniqueDBs(j)),1));
        end
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(i) & firstSpikeTimes(:,3) == uniqueDBs(j) & firstSpikeTimes(:,1) < spikeProbWindow(2) & firstSpikeTimes(:,1) > spikeProbWindow(1),1))
            firstSpikeStats(i,j,2) = size(find(firstSpikeTimes(:,2) == uniqueFreqs(i) &...
            firstSpikeTimes(:,3) == uniqueDBs(j) & firstSpikeTimes(:,1) < spikeProbWindow(2) & firstSpikeTimes(:,1) > spikeProbWindow(1)),1)/toneReps;
        end
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(i) & firstSpikeTimes(:,3) == uniqueDBs(j) & firstSpikeTimes(:,1) < spikeProbWindow(3) & firstSpikeTimes(:,1) > spikeProbWindow(2),1))
            firstSpikeStats(i,j,3) = size(find(firstSpikeTimes(:,2) == uniqueFreqs(i) &...
            firstSpikeTimes(:,3) == uniqueDBs(j) & firstSpikeTimes(:,1) < spikeProbWindow(3) & firstSpikeTimes(:,1) > spikeProbWindow(2)),1)/toneReps;
        end
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(i) & firstSpikeTimes(:,3) == uniqueDBs(j) & firstSpikeTimes(:,1) < spikeProbWindow(4) & firstSpikeTimes(:,1) > spikeProbWindow(3),1))
            firstSpikeStats(i,j,4) = size(find(firstSpikeTimes(:,2) == uniqueFreqs(i) &...
            firstSpikeTimes(:,3) == uniqueDBs(j) & firstSpikeTimes(:,1) < spikeProbWindow(4) & firstSpikeTimes(:,1) > spikeProbWindow(3)),1)/toneReps;
        end
        
    end
end

%% looking at spikes in general. this way i can look at multiple windows?

toneSpikeTimes = zeros(totalTrialNum,5);

%pulls out all the spikes per trial within defined bins.
for i = 1:totalTrialNum
    if ~isempty(find(rasterData(:,1) == i & rasterData(:,2) >0& rasterData(:,2) < spikeProbWindow(4),1))
        findToneSpikes = find(rasterData(:,1) == i & rasterData(:,2) >0 & rasterData(:,2) < spikeProbWindow(4));
        toneSpikes = rasterData(findToneSpikes,2);
        toneSpikeTimes(i,1) = size(toneSpikes(toneSpikes<spikeProbWindow(2)),1);
        toneSpikeTimes(i,2) = size(toneSpikes(toneSpikes<spikeProbWindow(3)&toneSpikes>spikeProbWindow(2)),1);
        toneSpikeTimes(i,3) = size(toneSpikes(toneSpikes<spikeProbWindow(4)&toneSpikes>spikeProbWindow(3)),1);
        toneSpikeTimes(i,4) = rasterData(findToneSpikes(1),3);
        toneSpikeTimes(i,5) = rasterData(findToneSpikes(1),4);
    end
end

%clears empty stuff
toneSpikeTimes(toneSpikeTimes(:,5) == 0,:) = [];

toneSpikeStats = zeros(length(uniqueFreqs),length(uniqueDBs),3);
for i =1:size(uniqueFreqs,1)
    for j = 1:size(uniqueDBs,1)
        if ~isempty(find(toneSpikeTimes(:,4) == uniqueFreqs(i)&toneSpikeTimes(:,5) == uniqueDBs(j),1))
            toneSpikeStats(i,j,1) = size(find(toneSpikeTimes(:,4) == uniqueFreqs(i) & toneSpikeTimes(:,5) == uniqueDBs(j) & toneSpikeTimes(:,1) > 0),1)/toneReps;
            toneSpikeStats(i,j,2) = size(find(toneSpikeTimes(:,4) == uniqueFreqs(i) & toneSpikeTimes(:,5) == uniqueDBs(j) & toneSpikeTimes(:,2) > 0),1)/toneReps;
            toneSpikeStats(i,j,3) = size(find(toneSpikeTimes(:,4) == uniqueFreqs(i) & toneSpikeTimes(:,5) == uniqueDBs(j) & toneSpikeTimes(:,3) > 0),1)/toneReps;
        end
    end
end


%alternate approach: calculate average firing rate based on entire 100ms
%bin before tone, find stats between trials. Tested, doesnt seem to work
%terribly well: huge variation! Q_Q
% binnedBaseline = zeros(1,totalTrialNum);
% for i = 1:totalTrialNum
%     binnedBaseline(i) = size(find(rasterData(:,1) == i & rasterData(:,2)<=0),1);
%     binnedBaseline(i) = binnedBaseline(i)*(preToneDuration(2)-preToneDuration(1));
% end