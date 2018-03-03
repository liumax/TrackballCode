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
    jumpsBack = round((s.Parameters.RasterWindow(1))/s.Parameters.InterpolationStepRotary);
    jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);

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
    
    hybridFind = sort([fullStopZeroFind,fullVelFind]);
    
    velFullStop(:,i) = mean(velRaster(:,fullStopZeroFind)');
    velFullVel(:,i) = mean(velRaster(:,fullVelFind)');
    velHybrid(:,i) = mean(velRaster(:,hybridFind)');
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
    trialsStore{i,5} = hybridFind;
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
        meansHolder(subsubInd,5) = nanmean((holder(hybridFind,3) - holder(hybridFind,1))./(holder(hybridFind,3) + holder(hybridFind,1)));
        meansHolder(subsubInd,6) = nanmean((holder(:,3)-holder(:,1))./(holder(:,3)+holder(:,1)));
        numHolder(subsubInd,1) = length(find(~isnan(((holder(fullVelFind,3) - holder(fullVelFind,1))./(holder(fullVelFind,3) + holder(fullVelFind,1))))));
        numHolder(subsubInd,2) = length(find(~isnan((holder(fullStopZeroFind,3) - holder(fullStopZeroFind,1))./(holder(fullStopZeroFind,3) + holder(fullStopZeroFind,1)))));
        numHolder(subsubInd,3) = length(find(~isnan((holder(laserZeroFind,3) - holder(laserZeroFind,1))./(holder(laserZeroFind,3) + holder(laserZeroFind,1)))));
        numHolder(subsubInd,4) = length(find(~isnan((holder(preZeroFind,3) - holder(preZeroFind,1))./(holder(preZeroFind,3) + holder(preZeroFind,1)))));
        numHolder(subsubInd,5) = length(find(~isnan((holder(hybridFind,3) - holder(hybridFind,1))./(holder(hybridFind,3) + holder(hybridFind,1)))));
        numHolder(subsubInd,6) = length(find(~isnan((holder(:,3)-holder(:,1))./(holder(:,3)+holder(:,1)))));
        numHolder(subsubInd,7) = numLaser;
        varHold(counter + subsubInd - 1) = var(holder(:,3));
        preMeanHold(counter+subsubInd-1) = mean(holder(:,3));
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
        histHybrid(counter+subsubInd-1,:) = mean(histStore(hybridFind,:))/histBin;
        histFull(counter+subsubInd-1,:) = mean(histStore);
        
        %calculate z score off all trials
        holder = (histStore)./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
%         baseHolder = reshape(baseHolder,1,[]); %180207 found this was
%         dramatically overinflating standard deviation. Since I'm only
%         looking at mean values, this makes no sense.
        baseHolder = mean(baseHolder);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        zhistFullVel(counter+subsubInd-1,:) = ((mean(holder(fullVelFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistFullStop(counter+subsubInd-1,:) = ((mean(holder(fullStopZeroFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistLaserStop(counter+subsubInd-1,:) = ((mean(holder(laserZeroFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistPreStop(counter+subsubInd-1,:) = ((mean(holder(preZeroFind,:))/histBin) - baseMean)/baseSTD;
        
        zhistFull(counter+subsubInd-1,:) = (mean(holder)-baseMean)/baseSTD;
        
        zhistHybrid(counter+subsubInd-1,:) = ((mean(holder(hybridFind,:))/histBin) - baseMean)/baseSTD;
        %calculate z scores individually
        holder = (histStore(fullVelFind,:))./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
%         baseHolder = reshape(baseHolder,1,[]);
        baseHolder = mean(baseHolder);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        speczhistFullVel(counter+subsubInd-1,:) = ((mean(histStore(fullVelFind,:))/histBin) - baseMean)/baseSTD;
        holder = (histStore(fullStopZeroFind,:))./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
%         baseHolder = reshape(baseHolder,1,[]);
        baseHolder = mean(baseHolder);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        speczhistFullStop(counter+subsubInd-1,:) = ((mean(histStore(fullStopZeroFind,:))/histBin) - baseMean)/baseSTD;


        holder = (histStore(hybridFind,:))./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
 %         baseHolder = reshape(baseHolder,1,[]);
        baseHolder = mean(baseHolder);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        speczhistHybrid(counter+subsubInd-1,:) = ((mean(histStore(hybridFind,:))/histBin) - baseMean)/baseSTD;

        holder = (histStore(laserZeroFind,:))./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
%         baseHolder = reshape(baseHolder,1,[]);
        baseHolder = mean(baseHolder);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        speczhistLaserStop(counter+subsubInd-1,:) = ((mean(histStore(laserZeroFind,:))/histBin) - baseMean)/baseSTD;
        holder = (histStore(preZeroFind,:))./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
%         baseHolder = reshape(baseHolder,1,[]);
        baseHolder = mean(baseHolder);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        speczhistPreStop(counter+subsubInd-1,:) = ((mean(histStore(preZeroFind,:))/histBin) - baseMean)/baseSTD;
        holder = (histStore(:,:))./histBin;
        zeroFind = find(histVector < 0,1,'last');
        baseHolder = holder(:,[1:zeroFind]);
%         baseHolder = reshape(baseHolder,1,[]);
        baseHolder = mean(baseHolder);
        baseMean = mean(baseHolder);
        baseSTD = std(baseHolder);
        speczhistFull(counter+subsubInd-1,:) = mean(histStore);
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

figure
plot(zhistFull')

%generate vector for velocity
velVector = [s.Parameters.RasterWindow(1):s.Parameters.InterpolationStepRotary:s.Parameters.RasterWindow(2)];
%identify msns and pvs
pvs = find(spikeWidthStore(:,1) < 4*10^-4);
msns = find(spikeWidthStore(:,1) >= 5*10^-4);

histVectMod = [-1:0.1:1];


%lets show what happens after the laser!
msnHists = histFull(msns,:);
preLaserWind = mean(msnHists(:,1:5)');
postLaserWind = mean(msnHists(:,40:45)');
msnzHists = zhistFull(msns,:);
preLaserWindZ = mean(msnzHists(:,1:5)');
postLaserWindZ = mean(msnzHists(:,40:45)');

hFig = figure
subplot(2,1,1)
hold on
for i = 1:length(msns)
    plot([preLaserWind(i),postLaserWind(i)],'c')
end
plot([mean(preLaserWind) mean(postLaserWind)],'k','LineWidth',2)
xlim([0 3])
xlabel('Pre vs Post Laser 500ms')
ylabel('Hz')
title('mean500ms bins at beginning of Pre and immediately after laser')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(2,1,2)
hold on
for i = 1:length(msns)
    plot([preLaserWindZ(i),postLaserWindZ(i)],'c')
end
plot([mean(preLaserWindZ) mean(postLaserWindZ)],'k','LineWidth',2)
xlim([0 3])
xlabel('Pre vs Post Laser 500ms')
ylabel('zScore')
title('mean500ms bins at beginning of Pre and immediately after laser')
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'preAndPostBin500ms';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%% Plot out histograms of modulation indices
%eliminate nan values
modMasterMeansHolder = masterMeansHolder;
modMasterIndex = [1:size(modMasterMeansHolder,1)];
for i = 1:size(modMasterMeansHolder,2)
    nanFind = find(isnan(modMasterMeansHolder(:,i)));
    modMasterMeansHolder(nanFind,:) = [];
    modMasterIndex(nanFind) = [];
end

[C ia ib] = intersect(modMasterIndex,msns);
modMasterMSNs = ia;
[C ia ib] = intersect(modMasterIndex,pvs);
modMasterPVs = ia;

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(4,1,1)
hist(modMasterMeansHolder(modMasterMSNs,1),histVectMod)
xlim([-1 1])
title('MSN Locomotor Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,2)
hist(modMasterMeansHolder(modMasterMSNs,2),histVectMod)
xlim([-1 1])
title('MSN Stationary Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,3)
hist(modMasterMeansHolder(modMasterMSNs,3),histVectMod)
xlim([-1 1])
title('MSN Laser Stop Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,4)
hist(modMasterMeansHolder(modMasterMSNs,4),histVectMod)
xlim([-1 1])
title('MSN Laser Start Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'LocomotionEffectsOnModIndex';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(4,1,1)
hist(modMasterMeansHolder(modMasterPVs,1),histVectMod)
xlim([-1 1])
title('PV Locomotor Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,2)
hist(modMasterMeansHolder(modMasterPVs,2),histVectMod)
xlim([-1 1])
title('PV Stationary Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,3)
hist(modMasterMeansHolder(modMasterPVs,3),histVectMod)
xlim([-1 1])
title('PV Laser Stop Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,4)
hist(modMasterMeansHolder(modMasterPVs,4),histVectMod)
xlim([-1 1])
title('PV Laser Start Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'LocomotionEffectsOnModIndexPV';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

% now lets plot out mean values as bar graphs?



hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
hold on
for i = 1:length(modMasterMSNs)
    plot(modMasterMeansHolder(modMasterMSNs(i),1:4),'k')
end
calcMean = nanmean(modMasterMeansHolder(modMasterMSNs,1:4));
calcSTE = nanstd(modMasterMeansHolder(modMasterMSNs,1:4))/sqrt(length(modMasterMSNs));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
% meanSig1 = signrank(modMasterMeansHolder(modMasterMSNs,1),modMasterMeansHolder(modMasterMSNs,4));
% meanSig2 = signrank(modMasterMeansHolder(modMasterMSNs,2),modMasterMeansHolder(modMasterMSNs,5));
% meanSig3 = signrank(modMasterMeansHolder(modMasterMSNs,3),modMasterMeansHolder(modMasterMSNs,5));
% title(strcat('Average Firing Rate modMasterMSNs (Hz) SIG34',num2str(meanSig1),'SIG45',num2str(meanSig2),'SIG35',num2str(meanSig3)))
plot([0.5 4.5],[0 0],'b')
xlim([0 5])
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(2,1,2)
hold on
for i = 1:length(modMasterPVs)
    plot(modMasterMeansHolder(modMasterPVs(i),1:4),'k')
end
calcMean = nanmean(modMasterMeansHolder(modMasterPVs,1:4));
calcSTE = nanstd(modMasterMeansHolder(modMasterPVs,1:4))/sqrt(length(modMasterPVs));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
% meanSig1 = signrank(modMasterMeansHolder(modMasterPVs,3),modMasterMeansHolder(modMasterPVs,4));
% meanSig2 = signrank(modMasterMeansHolder(modMasterPVs,4),modMasterMeansHolder(modMasterPVs,5));
% meanSig3 = signrank(modMasterMeansHolder(modMasterPVs,3),modMasterMeansHolder(modMasterPVs,5));
% title(strcat('Average Firing Rate modMasterPVs (Hz) SIG34',num2str(meanSig1),'SIG45',num2str(meanSig2),'SIG35',num2str(meanSig3)))
plot([0.5 4.5],[0 0],'b')
xlim([0 5])
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnANDpvVsStartVsStop';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets make the above with bar plots with error bars. 


hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
hold on
calcMean = nanmean(modMasterMeansHolder(modMasterMSNs,1:4));
calcSTE = nanstd(modMasterMeansHolder(modMasterMSNs,1:4))/sqrt(length(modMasterMSNs));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
bar(calcMean)
% meanSig1 = signrank(modMasterMeansHolder(modMasterMSNs,1),modMasterMeansHolder(modMasterMSNs,4));
% meanSig2 = signrank(modMasterMeansHolder(modMasterMSNs,2),modMasterMeansHolder(modMasterMSNs,5));
% meanSig3 = signrank(modMasterMeansHolder(modMasterMSNs,3),modMasterMeansHolder(modMasterMSNs,5));
% title(strcat('Average Firing Rate modMasterMSNs (Hz) SIG34',num2str(meanSig1),'SIG45',num2str(meanSig2),'SIG35',num2str(meanSig3)))
plot([0.5 4.5],[0 0],'b')
xlim([0.5 4.5])
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(2,1,2)
hold on
calcMean = nanmean(modMasterMeansHolder(modMasterPVs,1:4));
calcSTE = nanstd(modMasterMeansHolder(modMasterPVs,1:4))/sqrt(length(modMasterPVs));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
bar(calcMean)
% meanSig1 = signrank(modMasterMeansHolder(modMasterPVs,3),modMasterMeansHolder(modMasterPVs,4));
% meanSig2 = signrank(modMasterMeansHolder(modMasterPVs,4),modMasterMeansHolder(modMasterPVs,5));
% meanSig3 = signrank(modMasterMeansHolder(modMasterPVs,3),modMasterMeansHolder(modMasterPVs,5));
% title(strcat('Average Firing Rate modMasterPVs (Hz) SIG34',num2str(meanSig1),'SIG45',num2str(meanSig2),'SIG35',num2str(meanSig3)))
plot([0.5 4.5],[0 0],'b')
xlim([0.5 4.5])
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnANDpvVsStartVsStopBarPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%% now look at distribution of trial types
tester = typeStore;
for i = 1:4
tester(:,i) = tester(:,i)./tester(:,5);
end
%make pie chart of times in different states
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
% labels={'Loco','Stat','Lstop','Lgo'};
% pie(mean(tester(:,1:4)),labels)
pie(mean(tester(:,1:4)))
spikeGraphName = 'pieChartLocoState';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%also display as histograms
max(masterNumHolder(:,1:4));

hFig = figure

subplot(4,1,1)
hist(masterNumHolder(:,1),100)
% xlim([-1 1])
title('Running Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,2)
hist(masterNumHolder(:,2),100)
% xlim([-1 1])
title('Stationary Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,3)
hist(masterNumHolder(:,3),100)
% xlim([-1 1])
title('GotoStop Trials')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,4)
hist(masterNumHolder(:,4),100)
% xlim([-1 1])
title('StoptoGo Trials')
set(gca,'TickDir','out'); % The only other option is 'in'


%% lets eliminate all things with less than X trials, see what this changes

newMeansHolder = masterMeansHolder;

minTrialNum = 10;

shitFinder = find(masterNumHolder(:,1) <minTrialNum);
newMeansHolder(shitFinder,1) = NaN;

shitFinder = find(masterNumHolder(:,2) <minTrialNum);
newMeansHolder(shitFinder,2) = NaN;

shitFinder = find(masterNumHolder(:,3) <minTrialNum);
newMeansHolder(shitFinder,3) = NaN;

shitFinder = find(masterNumHolder(:,4) <minTrialNum);
newMeansHolder(shitFinder,4) = NaN;


hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(4,1,1)
hist(newMeansHolder(msns,1),histVectMod)
xlim([-1 1])
title('msns Running Trials elim <10')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,2)
hist(newMeansHolder(msns,2),histVectMod)
xlim([-1 1])
title('msns Stationary Trials elim <10')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,3)
hist(newMeansHolder(msns,3),histVectMod)
xlim([-1 1])
title('msns GotoStop Trials elim <10')
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(4,1,4)
hist(newMeansHolder(msns,4),histVectMod)
xlim([-1 1])
title('msns StoptoGo Trials elim <10')
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'modulationIndicesMSNelim10';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets re-index msns and pvs. 

pruneMeans = newMeansHolder(:,1:2);
pruneIndex = [1:length(pruneMeans)];
shitFinder = find(isnan(pruneMeans(:,1)));
pruneIndex(shitFinder) = [];
pruneMeans(shitFinder,:) = [];
shitFinder = find(isnan(pruneMeans(:,2)));
pruneIndex(shitFinder) = [];
pruneMeans(shitFinder,:) = [];


[C ia ib] = intersect(pruneIndex,msns);
pruneMSNs = ia;
[C ia ib] = intersect(pruneIndex,pvs);
prunePVs = ia;

%lets try this again, but with selection for only locomotion and stationary
%trials

simpMeansHolder = masterMeansHolder;

shitFinder = find(isnan(masterNumHolder(:,1)));
simpMeansHolder(shitFinder,1) = NaN;

shitFinder = find(isnan(masterNumHolder(:,2)));
simpMeansHolder(shitFinder,2) = NaN;


pruneMeans = simpMeansHolder(:,1:2);
pruneIndex = [1:length(pruneMeans)];
shitFinder = find(isnan(pruneMeans(:,1)));
pruneIndex(shitFinder) = [];
pruneMeans(shitFinder,:) = [];
shitFinder = find(isnan(pruneMeans(:,2)));
pruneIndex(shitFinder) = [];
pruneMeans(shitFinder,:) = [];


[C ia ib] = intersect(pruneIndex,msns);
pruneMSNs = ia;
[C ia ib] = intersect(pruneIndex,pvs);
prunePVs = ia;






%% plot out firing rate PSTHs for MSNs

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
tempHold = [];
for i = 1:length(spikeWidthStore)
    tempHold(:,i) = velFullStop(:,spikeWidthStore(i,2));
end
steHistFullStop = nanstd(tempHold');

adjHolder = velFullVel;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistFullVel = nansum(adjHolder,2);
tempHold = [];
for i = 1:length(spikeWidthStore)
    tempHold(:,i) = velFullVel(:,spikeWidthStore(i,2));
end
steHistFullVel = nanstd(tempHold');


adjHolder = velLaserStop;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistLaserStop = nansum(adjHolder,2);
tempHold = [];
for i = 1:length(spikeWidthStore)
    tempHold(:,i) = velLaserStop(:,spikeWidthStore(i,2));
end
steHistLaserStop = nanstd(tempHold');


adjHolder = velPreStop;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistPreStop = nansum(adjHolder,2);
tempHold = [];
for i = 1:length(spikeWidthStore)
    tempHold(:,i) = velPreStop(:,spikeWidthStore(i,2));
end
steHistPreStop = nanstd(tempHold');


adjHolder = velHybrid;
for i = 1:length(ratioStore)
    adjHolder(:,i) = adjHolder(:,i)*ratioStore(i);
end
meanHistHybrid = nansum(adjHolder,2);
tempHold = [];
for i = 1:length(spikeWidthStore)
    tempHold(:,i) = velHybrid(:,spikeWidthStore(i,2));
end
steHistHybrid = nanstd(tempHold');


hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hold on
plot(histVector,nanmean(histFullStop(msns,:)),'b-')
plot(histVector,nanmean(histFullStop(msns,:))+nanstd(histFullStop(msns,:))/sqrt(length(msns)),'b-')
plot(histVector,nanmean(histFullStop(msns,:))-nanstd(histFullStop(msns,:))/sqrt(length(msns)),'b-')

plot(histVector,nanmean(histFullVel(msns,:)),'r-')
plot(histVector,nanmean(histFullVel(msns,:))+nanstd(histFullVel(msns,:))/sqrt(length(msns)),'r-')
plot(histVector,nanmean(histFullVel(msns,:))-nanstd(histFullVel(msns,:))/sqrt(length(msns)),'r-')
title('Average Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'

%plot out Speeds
subplot(2,2,2)
hold on
plot(velVector,(meanHistFullStop),'b-')
plot(velVector,(meanHistFullStop)+steHistFullStop','b-')
plot(velVector,(meanHistFullStop)-steHistFullStop','b-')

plot(velVector,(meanHistFullVel),'r-')
plot(velVector,(meanHistFullVel)+steHistFullVel','r-')
plot(velVector,(meanHistFullVel)-steHistFullVel','r-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

subplot(2,2,3)
hold on
plot(histVector,nanmean(histLaserStop(msns,:)),'g-')
plot(histVector,nanmean(histLaserStop(msns,:))+nanstd(histLaserStop(msns,:))/sqrt(length(msns)),'g-')
plot(histVector,nanmean(histLaserStop(msns,:))-nanstd(histLaserStop(msns,:))/sqrt(length(msns)),'g-')
plot(histVector,nanmean(histPreStop(msns,:)),'c-')
plot(histVector,nanmean(histPreStop(msns,:))+nanstd(histPreStop(msns,:))/sqrt(length(msns)),'c-')
plot(histVector,nanmean(histPreStop(msns,:))-nanstd(histPreStop(msns,:))/sqrt(length(msns)),'c-')
title('Average Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'

%plot out Speeds
subplot(2,2,4)
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistLaserStop)+steHistLaserStop','g-')
plot(velVector,(meanHistLaserStop)-steHistLaserStop','g-')
plot(velVector,(meanHistPreStop),'c-')
plot(velVector,(meanHistPreStop)+steHistPreStop','c-')
plot(velVector,(meanHistPreStop)-steHistPreStop','c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnAverageRateVsVelocity';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%plot out firing rate stuffs by z score
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hold on
plot(histVector,nanmean(zhistFullStop(msns,:)),'b-')
plot(histVector,nanmean(zhistFullStop(msns,:))+nanstd(zhistFullStop(msns,:))/sqrt(length(msns)),'b-')
plot(histVector,nanmean(zhistFullStop(msns,:))-nanstd(zhistFullStop(msns,:))/sqrt(length(msns)),'b-')

plot(histVector,nanmean(zhistFullVel(msns,:)),'r-')
plot(histVector,nanmean(zhistFullVel(msns,:))+nanstd(zhistFullVel(msns,:))/sqrt(length(msns)),'r-')
plot(histVector,nanmean(zhistFullVel(msns,:))-nanstd(zhistFullVel(msns,:))/sqrt(length(msns)),'r-')

title('Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'

%plot out Speeds
subplot(2,2,2)
hold on
plot(velVector,(meanHistFullStop),'b-')
plot(velVector,(meanHistFullStop)+steHistFullStop','b-')
plot(velVector,(meanHistFullStop)-steHistFullStop','b-')

plot(velVector,(meanHistFullVel),'r-')
plot(velVector,(meanHistFullVel)+steHistFullVel','r-')
plot(velVector,(meanHistFullVel)-steHistFullVel','r-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

subplot(2,2,3)
hold on
plot(histVector,nanmean(zhistLaserStop(msns,:)),'g-')
plot(histVector,nanmean(zhistLaserStop(msns,:))+nanstd(zhistLaserStop(msns,:))/sqrt(length(msns)),'g-')
plot(histVector,nanmean(zhistLaserStop(msns,:))-nanstd(zhistLaserStop(msns,:))/sqrt(length(msns)),'g-')
plot(histVector,nanmean(zhistPreStop(msns,:)),'c-')
plot(histVector,nanmean(zhistPreStop(msns,:))+nanstd(zhistPreStop(msns,:))/sqrt(length(msns)),'c-')
plot(histVector,nanmean(zhistPreStop(msns,:))-nanstd(zhistPreStop(msns,:))/sqrt(length(msns)),'c-')
title('Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,4)
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistLaserStop)+steHistLaserStop','g-')
plot(velVector,(meanHistLaserStop)-steHistLaserStop','g-')
plot(velVector,(meanHistPreStop),'c-')
plot(velVector,(meanHistPreStop)+steHistPreStop','c-')
plot(velVector,(meanHistPreStop)-steHistPreStop','c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnAverageZVsVelocity';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now plot the z scores calculated of each individual trial
findInf = find(speczhistFullStop == Inf);
speczhistFullStop(findInf) = NaN;

findInf = find(speczhistFullVel == Inf);
speczhistFullVel(findInf) = NaN;

findInf = find(speczhistLaserStop == Inf);
speczhistLaserStop(findInf) = NaN;
findInf = find(speczhistLaserStop == -Inf);
speczhistLaserStop(findInf) = NaN;

findInf = find(speczhistPreStop == Inf);
speczhistPreStop(findInf) = NaN;
findInf = find(speczhistPreStop == -Inf);
speczhistPreStop(findInf) = NaN;

%plot out firing rate stuffs
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hold on
plot(histVector,nanmean(speczhistFullStop(msns,:)),'b-')
plot(histVector,nanmean(speczhistFullStop(msns,:))+nanstd(speczhistFullStop(msns,:))/sqrt(length(msns)),'b-')
plot(histVector,nanmean(speczhistFullStop(msns,:))-nanstd(speczhistFullStop(msns,:))/sqrt(length(msns)),'b-')
plot(histVector,nanmean(speczhistFullVel(msns,:)),'r-')
plot(histVector,nanmean(speczhistFullVel(msns,:))+nanstd(speczhistFullVel(msns,:))/sqrt(length(msns)),'r-')
plot(histVector,nanmean(speczhistFullVel(msns,:))-nanstd(speczhistFullVel(msns,:))/sqrt(length(msns)),'r-')
title('Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,2)
plot(velVector,(meanHistFullStop),'b-')
hold on
plot(velVector,(meanHistFullVel),'r-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

subplot(2,2,3)
hold on
plot(histVector,nanmean(speczhistLaserStop(msns,:)),'g-')
plot(histVector,nanmean(speczhistLaserStop(msns,:))+nanstd(speczhistLaserStop(msns,:))/sqrt(length(msns)),'g-')
plot(histVector,nanmean(speczhistLaserStop(msns,:))-nanstd(speczhistLaserStop(msns,:))/sqrt(length(msns)),'g-')
plot(histVector,nanmean(speczhistPreStop(msns,:)),'c-')
plot(histVector,nanmean(speczhistPreStop(msns,:))+nanstd(speczhistPreStop(msns,:))/sqrt(length(msns)),'c-')
plot(histVector,nanmean(speczhistPreStop(msns,:))-nanstd(speczhistPreStop(msns,:))/sqrt(length(msns)),'c-')
title('Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,4)
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistPreStop),'c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnAverageNewZVsVelocity';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Repeat for PVs
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hold on
plot(histVector,nanmean(histFullStop(pvs,:)),'b-')
plot(histVector,nanmean(histFullStop(pvs,:))+nanstd(histFullStop(pvs,:))/sqrt(length(pvs)),'b-')
plot(histVector,nanmean(histFullStop(pvs,:))-nanstd(histFullStop(pvs,:))/sqrt(length(pvs)),'b-')

plot(histVector,nanmean(histFullVel(pvs,:)),'r-')
plot(histVector,nanmean(histFullVel(pvs,:))+nanstd(histFullVel(pvs,:))/sqrt(length(pvs)),'r-')
plot(histVector,nanmean(histFullVel(pvs,:))-nanstd(histFullVel(pvs,:))/sqrt(length(pvs)),'r-')
title('Average PV Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,2)
plot(velVector,(meanHistFullStop),'b-')
hold on
plot(velVector,(meanHistFullVel),'r-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

subplot(2,2,3)
hold on
plot(histVector,nanmean(histLaserStop(pvs,:)),'g-')
plot(histVector,nanmean(histLaserStop(pvs,:))+nanstd(histLaserStop(pvs,:))/sqrt(length(pvs)),'g-')
plot(histVector,nanmean(histLaserStop(pvs,:))-nanstd(histLaserStop(pvs,:))/sqrt(length(pvs)),'g-')
plot(histVector,nanmean(histPreStop(pvs,:)),'c-')
plot(histVector,nanmean(histPreStop(pvs,:))+nanstd(histPreStop(pvs,:))/sqrt(length(pvs)),'c-')
plot(histVector,nanmean(histPreStop(pvs,:))-nanstd(histPreStop(pvs,:))/sqrt(length(pvs)),'c-')
title('Average PV Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,4)
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistPreStop),'c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'pvAverageRateVsVelocity';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%plot out firing rate stuffs
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hold on
plot(histVector,nanmean(zhistFullStop(pvs,:)),'b-')
plot(histVector,nanmean(zhistFullStop(pvs,:))+nanstd(zhistFullStop(pvs,:))/sqrt(length(pvs)),'b-')
plot(histVector,nanmean(zhistFullStop(pvs,:))-nanstd(zhistFullStop(pvs,:))/sqrt(length(pvs)),'b-')
plot(histVector,nanmean(zhistFullVel(pvs,:)),'r-')
plot(histVector,nanmean(zhistFullVel(pvs,:))+nanstd(zhistFullVel(pvs,:))/sqrt(length(pvs)),'r-')
plot(histVector,nanmean(zhistFullVel(pvs,:))-nanstd(zhistFullVel(pvs,:))/sqrt(length(pvs)),'r-')
title('Average PV Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,2)
plot(velVector,(meanHistFullStop),'b-')
hold on
plot(velVector,(meanHistFullVel),'r-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

subplot(2,2,3)
hold on
plot(histVector,nanmean(zhistLaserStop(pvs,:)),'g-')
plot(histVector,nanmean(zhistLaserStop(pvs,:))+nanstd(zhistLaserStop(pvs,:))/sqrt(length(pvs)),'g-')
plot(histVector,nanmean(zhistLaserStop(pvs,:))-nanstd(zhistLaserStop(pvs,:))/sqrt(length(pvs)),'g-')
plot(histVector,nanmean(zhistPreStop(pvs,:)),'c-')
plot(histVector,nanmean(zhistPreStop(pvs,:))+nanstd(zhistPreStop(pvs,:))/sqrt(length(pvs)),'c-')
plot(histVector,nanmean(zhistPreStop(pvs,:))-nanstd(zhistPreStop(pvs,:))/sqrt(length(pvs)),'c-')
title('Average PV Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,4)
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistPreStop),'c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow) 

set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'pvZRateVsVelocity';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

meanFRstore(:,1) = nanmean(histFullStop(:,1:zeroFind)');
meanFRstore(:,2) = nanmean(histFullVel(:,1:zeroFind)');
meanFRstore(:,3) = nanmean(histLaserStop(:,1:zeroFind)');
meanFRstore(:,4) = nanmean(histPreStop(:,1:zeroFind)');

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hold on
plot(histVector,nanmean(speczhistFullStop(pvs,:)),'b-')
plot(histVector,nanmean(speczhistFullStop(pvs,:))+nanstd(speczhistFullStop(pvs,:))/sqrt(length(pvs)),'b-')
plot(histVector,nanmean(speczhistFullStop(pvs,:))-nanstd(speczhistFullStop(pvs,:))/sqrt(length(pvs)),'b-')
plot(histVector,nanmean(speczhistFullVel(pvs,:)),'r-')
plot(histVector,nanmean(speczhistFullVel(pvs,:))+nanstd(speczhistFullVel(pvs,:))/sqrt(length(pvs)),'r-')
plot(histVector,nanmean(speczhistFullVel(pvs,:))-nanstd(speczhistFullVel(pvs,:))/sqrt(length(pvs)),'r-')
title('Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,2)
plot(velVector,(meanHistFullStop),'b-')
hold on
plot(velVector,(meanHistFullVel),'r-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

subplot(2,2,3)
hold on
plot(histVector,nanmean(speczhistLaserStop(pvs,:)),'g-')
plot(histVector,nanmean(speczhistLaserStop(pvs,:))+nanstd(speczhistLaserStop(pvs,:))/sqrt(length(pvs)),'g-')
plot(histVector,nanmean(speczhistLaserStop(pvs,:))-nanstd(speczhistLaserStop(pvs,:))/sqrt(length(pvs)),'g-')
plot(histVector,nanmean(speczhistPreStop(pvs,:)),'c-')
plot(histVector,nanmean(speczhistPreStop(pvs,:))+nanstd(speczhistPreStop(pvs,:))/sqrt(length(pvs)),'c-')
plot(histVector,nanmean(speczhistPreStop(pvs,:))-nanstd(speczhistPreStop(pvs,:))/sqrt(length(pvs)),'c-')
title('Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,2,4)
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistPreStop),'c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)

set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'pvAverageNewZVsVelocity';
savefig(hFig,spikeGraphName); 
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%% now lets plot out, for MSNS and PVs, hybrid vs increasing vs decreasing


%plot out firing rate stuffs
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
plot(histVector,nanmean(speczhistHybrid(msns,:)),'b-')
hold on
plot(histVector,nanmean(speczhistLaserStop(msns,:)),'g-')
plot(histVector,nanmean(speczhistPreStop(msns,:)),'c-')
title('MSN Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,1,2)
plot(velVector,(meanHistHybrid),'b-')
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistPreStop),'c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnsAverageNewZHybridVsVelocity';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
plot(histVector,nanmean(speczhistHybrid(pvs,:)),'b-')
hold on
plot(histVector,nanmean(speczhistLaserStop(pvs,:)),'g-')
plot(histVector,nanmean(speczhistPreStop(pvs,:)),'c-')
title('PV Average Z Firing Rates')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
%plot out Speeds
subplot(2,1,2)
plot(velVector,(meanHistHybrid),'b-')
hold on
plot(velVector,(meanHistLaserStop),'g-')
plot(velVector,(meanHistPreStop),'c-')
title('Average Velocities')
xlim(s.Parameters.RasterWindow)
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'pvsAverageNewZHybridVsVelocity';
savefig(hFig,spikeGraphName); 


%% now plot out average rate relative to stationary/locomotion
hFig = figure
subplot(2,1,1)
set(hFig, 'Position', [10 80 1240 850])
hold on
for i = 1:length(msns)
    plot(meanFRstore(msns(i),1:2),'k')
end
calcMean = nanmean(meanFRstore(msns,1:2));
calcSTE = nanstd(meanFRstore(msns,1:2))/sqrt(length(msns));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
meanSig = signrank(meanFRstore(msns,1),meanFRstore(msns,2));
title(strcat('Average Firing Rate msns (Hz) SIG',num2str(meanSig)))
xlim([0 3])
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(2,1,2)
hold on
for i = 1:length(pvs)
    plot(meanFRstore(pvs(i),1:2),'k')
end
calcMean = nanmean(meanFRstore(pvs,1:2));
calcSTE = nanstd(meanFRstore(pvs,1:2))/sqrt(length(pvs));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
meanSig = signrank(meanFRstore(pvs,1),meanFRstore(pvs,2));
title(strcat('Average Firing Rate pvs (Hz) SIG',num2str(meanSig)))
xlim([0 3])
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnPVaverageSpeedLocoStationary';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%Plot with only means
hFig = figure
subplot(2,1,1)
set(hFig, 'Position', [10 80 1240 850])
hold on

calcMean = nanmean(meanFRstore(msns,1:2));
calcSTE = nanstd(meanFRstore(msns,1:2))/sqrt(length(msns));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
meanSig = signrank(meanFRstore(msns,1),meanFRstore(msns,2));
title(strcat('Average Firing Rate msns (Hz) SIG',num2str(meanSig)))
xlim([0 3])
set(gca,'TickDir','out'); % The only other option is 'in'
subplot(2,1,2)
hold on

calcMean = nanmean(meanFRstore(pvs,1:2));
calcSTE = nanstd(meanFRstore(pvs,1:2))/sqrt(length(pvs));
errorbar(calcMean,calcSTE,'r','LineWidth',2)
meanSig = signrank(meanFRstore(pvs,1),meanFRstore(pvs,2));
title(strcat('Average Firing Rate pvs (Hz) SIG',num2str(meanSig)))
xlim([0 3])
set(gca,'TickDir','out'); % The only other option is 'in'
spikeGraphName = 'msnPVJustMeanaverageSpeedLocoStationary';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%calculate this with normalized firing rates
normMSN = meanFRstore(msns,2)./meanFRstore(msns,1);
normPV = meanFRstore(pvs,2)./meanFRstore(pvs,1);
findInf = find(normMSN == Inf | isnan(normMSN));
normMSN(findInf) = [];
findInf = find(normPV == Inf | isnan(normPV));
normPV(findInf) = [];

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(1,2,1)
errorbar([1,mean(normMSN)],[0,std(normMSN)/sqrt(length(msns))])

subplot(1,2,2)
errorbar([1,mean(normPV)],[0,std(normPV)/sqrt(length(msns))])










