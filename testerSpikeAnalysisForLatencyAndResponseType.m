histBin = 0.005; %hist bin in seconds
preToneDuration = [-0.1 0];
preToneBinCenters = [preToneDuration(1)+histBin/2:histBin:preToneDuration(2)-histBin/2];

toneHistBins = [0.01 0.005 0.001];
postToneDuration = [0 0.200];
% postToneBinCenters() = [postToneDuration(1)+toneHistBins/2:toneHistBins:postToneDuration(2)-toneHistBins/2]


%ONE DATA SET
load('160802_ML160721A_R_3300_IDfastTuningFullTuningAnalysis.mat') %loads analysis
load('160802_ML160721A_R_3300_IDfastTuning.mat') %loads sound file
%pulls raster data out
%fairly good tuned cell. good firing rate increase, but not very time
%locked
% rasterData = matclustStruct.matclust_param_nt7.Rasters{1,2};
%responsive cell, but extremely sparse. Very few spikes. Individual spikes
%pull z score to 4.0
% rasterData = matclustStruct.matclust_param_nt6.Rasters{1};
%responsive cell with low baseline firing. Extremely time locked response. 
rasterData = matclustStruct.matclust_param_nt10.Rasters{1};

%SECOND DATA SET
load('160802_ML160721A_R_3300_IDfastTuningFullTuningAnalysis.mat')
load('160802_ML160721A_R_3300_IDfastTuning.mat') %loads sound file



%plots to check for appropriate response
figure
plot(rasterData(:,2),rasterData(:,1),'b.')
%pulls frequency data
uniqueFreqs = unique(rasterData(:,3));
uniqueDBs = unique(rasterData(:,4));
toneReps = soundData.ToneRepetitions;
totalTrialNum = toneReps*size(uniqueFreqs,1)*size(uniqueDBs,1);

%one alternative to finding basline firing rate: average across all trials,
%then find mean and std from individual 5ms bins.
baselineHolder = rasterData(rasterData(:,2)<0,:);
%make histogram of baseline holder to generate baseline from which we can
%calculate response properties. 
[counts centers] = hist(baselineHolder(:,2),preToneBinCenters);
baselineCounts = counts/histBin/totalTrialNum;
baselineMean = mean(baselineCounts);
baselineSTD = std(baselineCounts);


%make a holder for rasters. this will allow for future analysis
rasterHolder = cell(size(uniqueFreqs,1),size(uniqueDBs,1));
histHolder = cell(size(uniqueFreqs,1),size(uniqueDBs,1));
sigHolder = zeros(size(uniqueFreqs,1),size(uniqueDBs,1),size(toneHistBins,2),2);
for i = 1:size(uniqueFreqs,1)
    for j = 1:size(uniqueDBs,1)
        rasterHolder{i,j} = rasterData(rasterData(:,3) == uniqueFreqs(i) & rasterData(:,4) == uniqueDBs(j) & rasterData(:,2)>0,:);
        histHolder{i,j} = cell(3,1);
        for k = 1:size(toneHistBins,2)
            postToneBinCenters = [postToneDuration(1)+toneHistBins(k)/2:toneHistBins(k):postToneDuration(2)-toneHistBins(k)/2];
            [counts centers] = hist(rasterHolder{i,j}(:,2),postToneBinCenters);
            histHolder{i,j}{k} = (counts/toneHistBins(k)/totalTrialNum-baselineMean)/baselineSTD;
            if ~isempty(find(histHolder{i,j}{k} > 3))
                sigHolder(i,j,k,1) = find(histHolder{i,j}{k} > 3,1,'first');
                sigHolder(i,j,k,2) = size(find(histHolder{i,j}{k} > 3),2);
            end
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