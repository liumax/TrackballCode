load('161019_ML161004A_L17_3000_IDfastTuningFullTuningAnalysis.mat')
figure
plot(s.nt19cluster1.AllHistograms)

testData = s.nt19cluster1;
testData = s.nt19cluster2;
rasterData = testData.AllRasters;
plot(rasterData(:,1),rasterData(:,2),'k.')

%pull negative values (pretone)
negDatas = rasterData(rasterData(:,1) < 0,1);
%sort these from negative to positive
sortNeg = sort(negDatas);
%calculate ISIs
negDiff = diff(sortNeg);

numReps = 500;
numSpikes = length(negDatas);
simDataHolder = zeros(numReps,numSpikes);
simDataHolder(:,1) = sortNeg(1); %seed with first spike time

tic
for i = 1:numReps
    %do random permutation of ISIs
    randInd = randperm(numSpikes-1);
    %now add these times to simDataHolder;
    simDataHolder(i,2:end) = negDiff(randInd);
    
end
toc

simDataHolder = cumsum(simDataHolder,2);

%generate histograms
simHist = hist(simDataHolder',[-0.1:0.005:0]);

%generate histogram from bin values
figure
hist(reshape(simHist,[],1),100)

[baselineCounts baselineCenters] = hist(reshape(simHist,[],1),200);

%really should just use the real data instead of binning it though. Here, I
%normalize the values to the number of trials and bin size

%% % NEED REPLACEMENT HERE%%%%% hardcoded values for testing
baselineData = reshape(simHist,[],1)/0.005/630;

%now i want to try and use this to calculate significant responses in the
%tone period. I do this by calculating percentiles
percentileRange = [95 99 99.9 99.99];
valueRange = zeros(length(percentileRange),1);

for i = 1:length(percentileRange)
    valueRange(i) = prctile(baselineData,percentileRange(i));
end

genHistData = testData.AllHistograms;
genHistData(:,2:length(percentileRange)+1) = zeros;

for i = 1:length(percentileRange)
    sigFinder = find(genHistData(:,1)>valueRange(i));
    genHistData(sigFinder,i+1) = 1;
end


figure
plot(1:length(testData.AllHistograms),genHistData(:,1))
hold on
plot(find(genHistData(:,2) == 1),genHistData(genHistData(:,2) == 1,1),'b*')
plot(find(genHistData(:,3) == 1),genHistData(genHistData(:,3) == 1,1),'g*')
plot(find(genHistData(:,4) == 1),genHistData(genHistData(:,4) == 1,1),'r*')



