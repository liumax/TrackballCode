load('160909_ML160718B_R_2200_8kHzPulseWithTermPulsedLaserPulsedSoundAnalysis.mat')


%want to pull binned number of spikes following tone presentation
binTimes = zeros(10,2);
binTimes(:,1) = [0:1:9];
binTimes(:,2) = binTimes(:,1) + 0.2;

numTrials = 30;

names = fields(matclustStruct);
names(end) = [];

cellCounter = 1;

massStore = zeros(1,60,10);


for j = 1:6
    numClust = matclustStruct.(names{j}).Clusters;
    for k = 1:numClust
        rasterData = matclustStruct.(names{j}).BigRasters{k};
        for i = 1:numTrials*2
            trialFinder = find(rasterData(:,1) == i);
%             laserCheck = find(rasterData(trialFinder,4) == 1);
            if ~isempty(trialFinder)
                for m = 1:10
                    spikeFinder = find(rasterData(trialFinder,2) > binTimes(m,1) & rasterData(trialFinder,2) < binTimes(m,2));
                    massStore(cellCounter,i,m) = length(spikeFinder);
                end
%             elseif ~isempty(trialFinder) & isempty(laserCheck)
                
            else
                
            end
        end
        cellCounter = cellCounter + 1;
    end
end

nonLaser = massStore(:,[1:2:end],:);
Laser = massStore(:,[2:2:end],:);


nonLaserMean = squeeze(mean(nonLaser,2));
LaserMean = squeeze(mean(Laser,2));

nonLaserSTD = squeeze(std(nonLaser,1,2))/sqrt(30);
LaserSTD = squeeze(std(Laser,1,2))/sqrt(30);