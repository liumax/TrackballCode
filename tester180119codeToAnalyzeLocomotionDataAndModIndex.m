%This is code to do locomotion analysis of data, and to find modulation
%index of running and non-running trials


targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);

counter = 1;
masterMeansHolder = [];
masterNumHolder = [];

for i = 1:numFiles
    load(targetFiles{i})

    dioTimes = s.Timing.LaserTimes;
    jumpsBack = round((s.Parameters.RasterWindow(1))/.01);
    jumpsForward = round(s.Parameters.RasterWindow(2)/.01);

    velRaster = zeros(jumpsForward-jumpsBack+1,length(dioTimes));
    for subsubInd = 1:length(dioTimes)
        %find the time from the velocity trace closest to the actual stim time
        targetInd = find(s.RotaryData.Velocity(:,1)-dioTimes(subsubInd) > 0,1,'first');
        %pull appropriate velocity data
        velRaster(:,subsubInd) = s.RotaryData.Velocity([targetInd+jumpsBack:targetInd+jumpsForward],2);
    end
    
    averageVel = mean(velRaster,2);

    preMean = mean(velRaster(1:-jumpsBack,:));
    preZeroFind = find(preMean == 0);
%     figure
%     plot(velRaster(:,preZeroFind))

    laserMean = mean(velRaster(-jumpsBack:-jumpsBack*2,:));
    laserZeroFind = find(laserMean == 0);
%     figure
%     plot(velRaster(:,laserZeroFind))

    [C ia ib] = intersect(preZeroFind,laserZeroFind);
    preZeroFind(ia) = [];
    laserZeroFind(ib) = [];

    fullStopMean = mean(velRaster(1:-jumpsBack*2,:));
    fullStopZeroFind = find(fullStopMean == 0);
%     figure
%     plot(velRaster(:,fullStopZeroFind))

    %find "running" trials
    combFind = sort([preZeroFind,laserZeroFind,fullStopZeroFind]);
    fullVelFind = [1:length(dioTimes)];
    fullVelFind(combFind) = [];
    
    velFullStop(:,i) = mean(velRaster(:,fullStopZeroFind)');
    velFullVel(:,i) = mean(velRaster(:,fullVelFind)');
    velLaserStop(:,i) = mean(velRaster(:,laserZeroFind)');
    velPreStop(:,i) = mean(velRaster(:,preZeroFind)');

    masterZero = zeros(length(dioTimes),4);
    masterZero(fullVelFind,1) = 1;
    masterZero(fullStopZeroFind,2) = 1;
    masterZero(laserZeroFind,3) = 1;
    masterZero(preZeroFind,4) = 1;
    bigMasterZero{i} = masterZero;
    disp('Stored Number of Each Trial Type')
    
    typeStore(i,:)  = [length(fullVelFind) length(fullStopZeroFind) length(laserZeroFind) length(preZeroFind) length(dioTimes)];
    trialsStore{i,1} = fullVelFind;
    trialsStore{i,2} = fullStopZeroFind;
    trialsStore{i,3} = laserZeroFind;
    trialsStore{i,4} = preZeroFind;
    numLaser = length(dioTimes);

    numUnits = length(s.DesignationName);
    meansHolder = [];
    numHolder = [];
    for subsubInd = 1:numUnits
        holder = s.(s.DesignationName{subsubInd}).TrialBinnedSpikes;
        meansHolder(subsubInd,1) = nanmean((holder(fullVelFind,3) - holder(fullVelFind,1))./(holder(fullVelFind,3) + holder(fullVelFind,1)));
        meansHolder(subsubInd,2) = nanmean((holder(fullStopZeroFind,3) - holder(fullStopZeroFind,1))./(holder(fullStopZeroFind,3) + holder(fullStopZeroFind,1)));
        meansHolder(subsubInd,3) = nanmean((holder(laserZeroFind,3) - holder(laserZeroFind,1))./(holder(laserZeroFind,3) + holder(laserZeroFind,1)));
        meansHolder(subsubInd,4) = nanmean((holder(preZeroFind,3) - holder(preZeroFind,1))./(holder(preZeroFind,3) + holder(preZeroFind,1)));
        meansHolder(subsubInd,5) = nanmean((holder(:,3)-holder(:,1))./(holder(:,3)+holder(:,1)));
        numHolder(subsubInd,1) = length(find(~isnan(((holder(fullVelFind,3) - holder(fullVelFind,1))./(holder(fullVelFind,3) + holder(fullVelFind,1))))));
        numHolder(subsubInd,2) = length(find(~isnan((holder(fullStopZeroFind,3) - holder(fullStopZeroFind,1))./(holder(fullStopZeroFind,3) + holder(fullStopZeroFind,1)))));
        numHolder(subsubInd,3) = length(find(~isnan((holder(laserZeroFind,3) - holder(laserZeroFind,1))./(holder(laserZeroFind,3) + holder(laserZeroFind,1)))));
        numHolder(subsubInd,4) = length(find(~isnan((holder(preZeroFind,3) - holder(preZeroFind,1))./(holder(preZeroFind,3) + holder(preZeroFind,1)))));
        numHolder(subsubInd,5) = length(find(~isnan((holder(:,3)-holder(:,1))./(holder(:,3)+holder(:,1)))));
        numHolder(subsubInd,6) = numLaser;
        %make and save psths of firing rates as well
        histBin = 0.1;
        histVector = [s.Parameters.RasterWindow(1)+histBin/2:histBin:s.Parameters.RasterWindow(2)-+histBin/2];
        histStore = zeros(length(dioTimes),length(histVector));
        for sousInd = 1:length(dioTimes)
            %check for spikes
            spikeFinder = find(s.(s.DesignationName{subsubInd}).RasterLaser(:,2) == sousInd);
            if spikeFinder
                histStore(sousInd,:) = hist(s.(s.DesignationName{subsubInd}).RasterLaser(spikeFinder,1),histVector);
            end
        end
        %store average values
        histFullVel(counter+subsubInd-1,:) = mean(histStore(fullVelFind,:))/histBin;
        histFullStop(counter+subsubInd-1,:) = mean(histStore(fullStopZeroFind,:))/histBin;
        histLaserStop(counter+subsubInd-1,:) = mean(histStore(laserZeroFind,:))/histBin;
        histPreStop(counter+subsubInd-1,:) = mean(histStore(preZeroFind,:))/histBin;
        histFull(counter+subsubInd-1,:) = mean(histStore);
        holder = (histStore(:,:))./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
        baseHolder = reshape(baseHolder,1,[]);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        zhistFullVel(counter+subsubInd-1,:) = ((mean(histStore(fullVelFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistFullStop(counter+subsubInd-1,:) = ((mean(histStore(fullStopZeroFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistLaserStop(counter+subsubInd-1,:) = ((mean(histStore(laserZeroFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistPreStop(counter+subsubInd-1,:) = ((mean(histStore(preZeroFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistFull(counter+subsubInd-1,:) = mean(histStore);
    end
    
    %store pv/msn designations
    spikeWidthStore(counter:counter+numUnits-1,1) = s.MasterSheet(:,2);
    averageRateStore(counter:counter+numUnits-1,1) = s.MasterSheet(:,3);
    spikeWidthStore(counter:counter+numUnits-1,2) = i;
    
    disp('Storing Mod Indices')
    masterMeansHolder(counter:counter+numUnits-1,:) = meansHolder;
    masterNumHolder(counter:counter+numUnits-1,:) = numHolder;
    disp('Updating Counter')
    counter = counter + numUnits;

end
disp(counter)

pvs = find(spikeWidthStore(:,1) < 4*10^-4);
msns = find(spikeWidthStore(:,1) >= 5*10^-4);

histVectMod = [-1:0.1:1];

figure
subplot(4,1,1)
hist(masterMeansHolder(msns,1),histVectMod)
xlim([-1 1])
title('msns Running Trials')

subplot(4,1,2)
hist(masterMeansHolder(msns,2),histVectMod)
xlim([-1 1])
title('msns Stationary Trials')

subplot(4,1,3)
hist(masterMeansHolder(msns,3),histVectMod)
xlim([-1 1])
title('msns GotoStop Trials')

subplot(4,1,4)
hist(masterMeansHolder(msns,4),histVectMod)
xlim([-1 1])
title('msns StoptoGo Trials')

%for anatol

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(4,1,1)
hist(masterMeansHolder(msns,1),histVectMod)
xlim([-1 1])
title('MSN Locomotor Trials')

subplot(4,1,2)
hist(masterMeansHolder(msns,2),histVectMod)
xlim([-1 1])
title('MSN Stationary Trials')

subplot(4,1,3)
hist(masterMeansHolder(msns,3),histVectMod)
xlim([-1 1])
title('MSN Laser Stop Trials')

subplot(4,1,4)
hist(masterMeansHolder(msns,4),histVectMod)
xlim([-1 1])
title('MSN Laser Start Trials')

spikeGraphName = 'LocomotionEffectsOnModIndex';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



figure
subplot(4,1,1)
hist(masterMeansHolder(pvs,1),histVectMod)
xlim([-1 1])
title('pvs Running Trials')

subplot(4,1,2)
hist(masterMeansHolder(pvs,2),histVectMod)
xlim([-1 1])
title('pvs Stationary Trials')

subplot(4,1,3)
hist(masterMeansHolder(pvs,3),histVectMod)
xlim([-1 1])
title('pvs GotoStop Trials')

subplot(4,1,4)
hist(masterMeansHolder(pvs,4),histVectMod)
xlim([-1 1])
title('pvs StoptoGo Trials')

%now look at distribution of trial types
tester = typeStore;
for i = 1:4
tester(:,i) = tester(:,i)./tester(:,5);
end


max(masterNumHolder(:,1:4));

figure
subplot(4,1,1)
hist(masterNumHolder(:,1),100)
% xlim([-1 1])
title('Running Trials')

subplot(4,1,2)
hist(masterNumHolder(:,2),100)
% xlim([-1 1])
title('Stationary Trials')

subplot(4,1,3)
hist(masterNumHolder(:,3),100)
% xlim([-1 1])
title('GotoStop Trials')

subplot(4,1,4)
hist(masterNumHolder(:,4),100)
% xlim([-1 1])
title('StoptoGo Trials')

%lets eliminate all things with less than 10 trials

newMeansHolder = masterMeansHolder;

shitFinder = find(masterNumHolder(:,1) <10);
newMeansHolder(shitFinder,1) = NaN;

shitFinder = find(masterNumHolder(:,2) <10);
newMeansHolder(shitFinder,2) = NaN;

shitFinder = find(masterNumHolder(:,3) <10);
newMeansHolder(shitFinder,3) = NaN;

shitFinder = find(masterNumHolder(:,4) <10);
newMeansHolder(shitFinder,4) = NaN;


figure
subplot(4,1,1)
hist(newMeansHolder(msns,1),histVectMod)
xlim([-1 1])
title('msns Running Trials')

subplot(4,1,2)
hist(newMeansHolder(msns,2),histVectMod)
xlim([-1 1])
title('msns Stationary Trials')

subplot(4,1,3)
hist(newMeansHolder(msns,3),histVectMod)
xlim([-1 1])
title('msns GotoStop Trials')

subplot(4,1,4)
hist(newMeansHolder(msns,4),histVectMod)
xlim([-1 1])
title('msns StoptoGo Trials')


%plot out firing rate stuffs

%first, need to properly proportion velocity traces
for i = 1:numFiles
    lengthFind = length(find(spikeWidthStore(:,2) == i));
    ratioStore(i) = lengthFind/(counter - 1);
end

adjHolder = velFullStop;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistFullStop = nansum(adjHolder,2);

adjHolder = velFullVel;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistFullVel = nansum(adjHolder,2);

adjHolder = velLaserStop;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistLaserStop = nansum(adjHolder,2);

adjHolder = velPreStop;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistPreStop = nansum(adjHolder,2);


figure
subplot(2,1,1)
plot(nanmean(histFullStop(msns,:)),'b.-')
hold on
plot(nanmean(histFullVel(msns,:)),'r.-')
plot(nanmean(histLaserStop(msns,:)),'g.-')
plot(nanmean(histPreStop(msns,:)),'c.-')
title('Average Firing Rates')
xlim([0 length(nanmean(histFullStop))])

%plot out Speeds
subplot(2,1,2)
plot((meanHistFullStop),'b.-')
hold on
plot((meanHistFullVel),'r.-')
plot((meanHistLaserStop),'g.-')
plot((meanHistPreStop),'c.-')
title('Average Velocities')
xlim([0 length(nanmean(velFullStop'))])


%plot out firing rate stuffs
figure
subplot(2,1,1)
plot(nanmean(zhistFullStop(msns,:)),'b.-')
hold on
plot(nanmean(zhistFullVel(msns,:)),'r.-')
plot(nanmean(zhistLaserStop(msns,:)),'g.-')
plot(nanmean(zhistPreStop(msns,:)),'c.-')
title('Average Z Firing Rates')
xlim([0 length(nanmean(zhistFullStop))])

%plot out Speeds
subplot(2,1,2)
plot((meanHistFullStop),'b.-')
hold on
plot((meanHistFullVel),'r.-')
plot((meanHistLaserStop),'g.-')
plot((meanHistPreStop),'c.-')
title('Average Velocities')
xlim([0 length(nanmean(velFullStop'))])
