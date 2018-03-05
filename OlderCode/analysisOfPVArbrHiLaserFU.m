



% clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);

%extract distances
for i = 1:numFiles
    targetName = targetFiles{i};
    findSpace = strfind(targetName,'_');
    probeDepth(i) = str2num(targetName(findSpace(3)+1:findSpace(4)-1));
end

masterInd = 1;
for i = 1:numFiles
   %load the target file
    load(targetFiles{i})
    numUnits = length(s.DesignationName);
    master(masterInd:masterInd + numUnits-1,:) = s.MasterSheet;
    for j = 1:numUnits
        holder(masterInd + j -1,1) = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(101:200,4));
        holder(masterInd + j -1,2) = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(101:200,5));
        holder(masterInd + j -1,3) = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(101:200,6));
        holder(masterInd + j -1,4) = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(1:100,4));
        holder(masterInd + j -1,5) = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(1:100,5));
        holder(masterInd + j -1,6) = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(1:100,6));
        infoStore(masterInd + j - 1) = mean(s.(s.DesignationName{j}).HistogramLaser(1:20));
    end
    masterInd = masterInd + numUnits;
    
    
    
end

master(:,26:31) = holder;

mod1 = (master(:,31)-master(:,29))./(master(:,31)+master(:,29));
mod2 = (master(:,28)-master(:,26))./(master(:,28)+master(:,26));

rateChange1 = master(:,31)-master(:,29);
rateChange2 = master(:,28)-master(:,26);


pvs = find(master(:,1) ==1);
msns = find(master(:,1) == 0);

medians(1) = median(mod1(pvs));
medians(2) = median(mod1(msns));
medians(3) = median(mod2(pvs));
medians(4) = median(mod2(msns));

pVals(1) = signrank(mod1(pvs));
pVals(2) = signrank(mod1(msns));
pVals(3) = signrank(mod2(pvs));
pVals(4) = signrank(mod2(msns));

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hist(mod1(pvs),[-1:0.1:1])
title(strcat('PVMod3mW Median:',num2str(medians(1)),'pVal:',num2str(pVals(1))))
xlim([-1 1])
ylim([0 8])

subplot(2,2,2)
hist(mod1(msns),[-1:0.1:1])
title(strcat('MSNMod3mW Median:',num2str(medians(2)),'pVal:',num2str(pVals(2))))
xlim([-1 1])
ylim([0 9])

subplot(2,2,3)
hist(mod2(pvs),[-1:0.1:1])
title(strcat('PVMod15mW Median:',num2str(medians(3)),'pVal:',num2str(pVals(3))))
xlim([-1 1])
ylim([0 8])

subplot(2,2,4)
hist(mod2(msns),[-1:0.1:1])
title(strcat('MSNMod15mW Median:',num2str(medians(4)),'pVal:',num2str(pVals(4))))
xlim([-1 1])
ylim([0 9])


spikeGraphName = 'ModulationHistograms';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now make a scatter plot with lines

%do stats test on change in modulation index
pvalmsn = signrank(mod1(msns),mod2(msns));
pvalpv = signrank(mod1(pvs),mod2(pvs));


hFig= figure;
subplot(2,1,1)
hold on
for i = 1:length(msns)
    plot([mod1(msns(i)) mod2(msns(i))],'k')
end
plot(ones(length(msns),1),mod1(msns),'b.')
plot(ones(length(msns),1)*2,mod2(msns),'g.')
errorbar([mean(mod1(msns)) mean(mod2(msns))],[std(mod1(msns))/sqrt(length(msns)) std(mod2(msns))/sqrt(length(msns))],'LineWidth',2)
xlim([0 3])
ylim([-1 1])
xlabel('3 mW vs 15 mW')
ylabel('Modulation Index')
title(strcat('msns Mod Index. Pval:',num2str(pvalmsn)))

subplot(2,1,2)
hold on
for i = 1:length(pvs)
    plot([mod1(pvs(i)) mod2(pvs(i))],'k')
end
plot(ones(length(pvs),1),mod1(pvs),'b.')
plot(ones(length(pvs),1)*2,mod2(pvs),'g.')
errorbar([mean(mod1(pvs)) mean(mod2(pvs))],[std(mod1(pvs))/sqrt(length(pvs)) std(mod2(pvs))/sqrt(length(pvs))],'LineWidth',2)
xlim([0 3])
ylim([-1 1])
xlabel('3 mW vs 15 mW')
ylabel('Modulation Index')
title(strcat('msns Mod Index. Pval:',num2str(pvalpv)))

spikeGraphName = 'pvAndmsnsModCompare';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%now lets look at the rateChange vs laser power
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
hold on
plot(1,rateChange1(pvs),'k.')
plot(1,mean(rateChange1(pvs)),'r.')
plot(2,rateChange2(pvs),'k.')
plot(2,mean(rateChange2(pvs)),'r.')
xlim([0 3])
subplot(2,1,2)
hold on
plot(1,rateChange1(msns),'k.')
plot(1,mean(rateChange1(msns)),'r.')
plot(2,rateChange2(msns),'k.')
plot(2,mean(rateChange2(msns)),'r.')
xlim([0 3])

%now lets pull out the average histograms
masterInd = 1;
for i = 1:numFiles
   %load the target file
    load(targetFiles{i})
    numUnits = length(s.DesignationName);
%     normHistStore = zeros(120,1,2);
%     zHistStore = zeros(120,1,2);
    %pull locomotion data
%     tempLoco = zeros(200,1);
    locoStep = s.Parameters.InterpolationStepRotary;
    rastWind = s.Parameters.RasterWindow;
    locoWind = round(rastWind/locoStep);
    locoStore = [];
    for j = 1:length(s.Timing.LaserTimes)
        findNear = find(s.RotaryData.Velocity(:,1) - s.Timing.LaserTimes(j) > 0,1,'first');
        if findNear + locoWind(2) < length(s.RotaryData.Velocity)
            locoStore(:,j) = s.RotaryData.Velocity(findNear+locoWind(1):findNear+locoWind(2),2);
        else
            locoStore(:,j) = NaN(length([locoWind(1):locoWind(2)]),1);
        end
    end
    bigLoco(:,(200*(i-1))+1:i*200) = locoStore;
    locoLow(:,(100*(i-1))+1:i*100) = locoStore(:,1:100);
    locoHi(:,(100*(i-1))+1:i*100) = locoStore(:,101:200);
    for j = 1:numUnits
        normHistStore(:,masterInd + j - 1,1) = mean(s.(s.DesignationName{j}).TrialHists(:,1:100)');
        normHistStore(:,masterInd + j - 1,2) = mean(s.(s.DesignationName{j}).TrialHists(:,101:200)');
        baselineVal = mean(reshape(s.(s.DesignationName{j}).TrialHists(1:30,1:100),[],1));
        baselineSTD = std(reshape(s.(s.DesignationName{j}).TrialHists(1:30,1:100),[],1));
        zHistStore(:,masterInd+j-1,1)= (mean(s.(s.DesignationName{j}).TrialHists(:,1:100)')-baselineVal)/baselineSTD;
        baselineVal = mean(reshape(s.(s.DesignationName{j}).TrialHists(1:30,101:200),[],1));
        baselineSTD = std(reshape(s.(s.DesignationName{j}).TrialHists(1:30,101:200),[],1));
        zHistStore(:,masterInd+j-1,2)= (mean(s.(s.DesignationName{j}).TrialHists(:,101:200)')-baselineVal)/baselineSTD;
        locoCount(masterInd) = i;
    end
    masterInd = masterInd + numUnits;
    
    
    
end

%lets do some basic locomotion analysis

dioTimes = s.Timing.LaserTimes;
jumpsBack = round((s.Parameters.RasterWindow(1))/s.Parameters.InterpolationStepRotary);
jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);

averageVel = mean(locoHi,2);

preMean = mean(locoHi(1:-jumpsBack,:));
preZeroFind = find(preMean == 0);
%     figure
%     plot(locoHi(:,preZeroFind))

laserMean = mean(locoHi(-jumpsBack:-jumpsBack*2,:));
laserZeroFind = find(laserMean == 0);
%     figure
%     plot(locoHi(:,laserZeroFind))

[C ia ib] = intersect(preZeroFind,laserZeroFind);
preZeroFind(ia) = [];
laserZeroFind(ib) = [];

fullStopMean = mean(locoHi(1:-jumpsBack*2,:));
fullStopZeroFind = find(fullStopMean == 0);
%     figure
%     plot(locoHi(:,fullStopZeroFind))

%find "running" trials
combFind = sort([preZeroFind,laserZeroFind,fullStopZeroFind]);
fullVelFind = [1:600];
fullVelFind(combFind) = [];

hybridFind = sort([fullStopZeroFind,fullVelFind]);

velFullStop(:,i) = mean(locoHi(:,fullStopZeroFind)');
velFullVel(:,i) = mean(locoHi(:,fullVelFind)');
velHybrid(:,i) = mean(locoHi(:,hybridFind)');
velLaserStop(:,i) = mean(locoHi(:,laserZeroFind)');
velPreStop(:,i) = mean(locoHi(:,preZeroFind)');

hFig = figure
subplot(2,1,1)

pie([length(fullVelFind),length(fullStopZeroFind),length(preZeroFind),length(laserZeroFind)],{'Loco','Stop','LaserGo','LaserStop'})
title('Hi Laser Power Trials')
averageVel = mean(locoLow,2);

preMean = mean(locoLow(1:-jumpsBack,:));
preZeroFind = find(preMean == 0);
%     figure
%     plot(locoLow(:,preZeroFind))

laserMean = mean(locoLow(-jumpsBack:-jumpsBack*2,:));
laserZeroFind = find(laserMean == 0);
%     figure
%     plot(locoLow(:,laserZeroFind))

[C ia ib] = intersect(preZeroFind,laserZeroFind);
preZeroFind(ia) = [];
laserZeroFind(ib) = [];

fullStopMean = mean(locoLow(1:-jumpsBack*2,:));
fullStopZeroFind = find(fullStopMean == 0);
%     figure
%     plot(locoLow(:,fullStopZeroFind))

%find "running" trials
combFind = sort([preZeroFind,laserZeroFind,fullStopZeroFind]);
fullVelFind = [1:600];
fullVelFind(combFind) = [];

hybridFind = sort([fullStopZeroFind,fullVelFind]);

velFullStop(:,i) = mean(locoLow(:,fullStopZeroFind)');
velFullVel(:,i) = mean(locoLow(:,fullVelFind)');
velHybrid(:,i) = mean(locoLow(:,hybridFind)');
velLaserStop(:,i) = mean(locoLow(:,laserZeroFind)');
velPreStop(:,i) = mean(locoLow(:,preZeroFind)');

subplot(2,1,2)

pie([length(fullVelFind),length(fullStopZeroFind),length(preZeroFind),length(laserZeroFind)],{'Loco','Stop','LaserGo','LaserStop'})
title('Low Laser Power Trials')

spikeGraphName = 'pieChartsLocoDist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')




rasterVect = [-2:0.05:4-0.05];

pvmin = min([mean(normHistStore(:,pvs,1)') mean(normHistStore(:,pvs,2)')])*20*.5;
pvmax = max([mean(normHistStore(:,pvs,1)') mean(normHistStore(:,pvs,2)')])*20*1.5;
msnmin = min([mean(normHistStore(:,msns,1)') mean(normHistStore(:,msns,2)')])*20*.5;
msnmax = max([mean(normHistStore(:,msns,1)') mean(normHistStore(:,msns,2)')])*20*1.5;

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
plot(rasterVect,mean(normHistStore(:,pvs,1)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(normHistStore(:,pvs,1)')*20 + std(normHistStore(:,pvs,1)'*20)/sqrt(length(pvs)));
plot(rasterVect,mean(normHistStore(:,pvs,1)')*20 - std(normHistStore(:,pvs,1)'*20)/sqrt(length(pvs)));
ylim([pvmin pvmax])
title('PVMod3mW Mean')

subplot(2,2,2)
% plot(rasterVect,mean(normHistStore(:,msns,1)')*20)
plot(rasterVect,mean(normHistStore(:,msns,1)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(normHistStore(:,msns,1)')*20 + std(normHistStore(:,msns,1)'*20)/sqrt(length(msns)));
plot(rasterVect,mean(normHistStore(:,msns,1)')*20 - std(normHistStore(:,msns,1)'*20)/sqrt(length(msns)));
title('MSNMod3mW Mean')
ylim([msnmin msnmax])

subplot(2,2,3)
plot(rasterVect,mean(normHistStore(:,pvs,2)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(normHistStore(:,pvs,2)')*20 + std(normHistStore(:,pvs,2)'*20)/sqrt(length(pvs)));
plot(rasterVect,mean(normHistStore(:,pvs,2)')*20 - std(normHistStore(:,pvs,2)'*20)/sqrt(length(pvs)));
% plot(rasterVect,mean(normHistStore(:,pvs,2)')*20)
title('PVMod15mW Mean')
ylim([pvmin pvmax])

subplot(2,2,4)
% plot(rasterVect,mean(normHistStore(:,msns,2)')*20)
plot(rasterVect,mean(normHistStore(:,msns,2)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(normHistStore(:,msns,2)')*20 + std(normHistStore(:,msns,2)'*20)/sqrt(length(msns)));
plot(rasterVect,mean(normHistStore(:,msns,2)')*20 - std(normHistStore(:,msns,2)'*20)/sqrt(length(msns)));
title('MSNMod15mW Mean')
ylim([msnmin msnmax])
spikeGraphName = 'Average Rates (Hz)';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



rasterVect = [-2:0.05:4-0.05];

pvmin = min([mean(zHistStore(:,pvs,1)') mean(zHistStore(:,pvs,2)')])*20*.5;
pvmax = max([mean(zHistStore(:,pvs,1)') mean(zHistStore(:,pvs,2)')])*20*1.5;
msnmin = min([mean(zHistStore(:,msns,1)') mean(zHistStore(:,msns,2)')])*20*.5;
msnmax = max([mean(zHistStore(:,msns,1)') mean(zHistStore(:,msns,2)')])*20*1.5;

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
plot(rasterVect,mean(zHistStore(:,pvs,1)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(zHistStore(:,pvs,1)')*20 + std(zHistStore(:,pvs,1)'*20)/sqrt(length(pvs)));
plot(rasterVect,mean(zHistStore(:,pvs,1)')*20 - std(zHistStore(:,pvs,1)'*20)/sqrt(length(pvs)));

title('PVMod3mW Mean')

subplot(2,2,2)
% plot(rasterVect,mean(normHistStore(:,msns,1)')*20)
plot(rasterVect,mean(zHistStore(:,msns,1)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(zHistStore(:,msns,1)')*20 + std(zHistStore(:,msns,1)'*20)/sqrt(length(msns)));
plot(rasterVect,mean(zHistStore(:,msns,1)')*20 - std(zHistStore(:,msns,1)'*20)/sqrt(length(msns)));
title('MSNMod3mW Mean')

subplot(2,2,3)
plot(rasterVect,mean(zHistStore(:,pvs,2)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(zHistStore(:,pvs,2)')*20 + std(zHistStore(:,pvs,2)'*20)/sqrt(length(pvs)));
plot(rasterVect,mean(zHistStore(:,pvs,2)')*20 - std(zHistStore(:,pvs,2)'*20)/sqrt(length(pvs)));
% plot(rasterVect,mean(normHistStore(:,pvs,2)')*20)
title('PVMod15mW Mean')

subplot(2,2,4)
% plot(rasterVect,mean(normHistStore(:,msns,2)')*20)
plot(rasterVect,mean(zHistStore(:,msns,2)')*20,'LineWidth',2)
hold on
plot(rasterVect,mean(zHistStore(:,msns,2)')*20 + std(zHistStore(:,msns,2)'*20)/sqrt(length(msns)));
plot(rasterVect,mean(zHistStore(:,msns,2)')*20 - std(zHistStore(:,msns,2)'*20)/sqrt(length(msns)));
title('MSNMod15mW Mean')

spikeGraphName = 'Average Rates (Z)';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at the early and late periods, based on z score
testPre = squeeze(mean(msnZHist(1:10,:,:)));
testPost = squeeze(mean(msnZHist(80:90,:,:)));

hFig = figure;
subplot(2,1,1)
hold on
for i = 1:length(msns)
    plot([testPre(i,1),testPost(i,1)],'c')
end
plot([mean(testPre(:,1)),mean(testPost(:,1))],'k','LineWidth',2)
xlim([0 3])

subplot(2,1,2)
hold on
for i = 1:length(msns)
    plot([testPre(i,2),testPost(i,2)],'c')
end
plot([mean(testPre(:,2)),mean(testPost(:,2))],'k','LineWidth',2)
xlim([0 3])

%try with spikes, not Z
testPre = squeeze(mean(normHistStore(1:10,msns,:)));
testPost = squeeze(mean(normHistStore(80:90,msns,:)));
modLow = (testPost(:,1) - testPre(:,1))./(testPost(:,1) + testPre(:,1));
modHi = (testPost(:,2) - testPre(:,2))./(testPost(:,2) + testPre(:,2));

hFig = figure;
subplot(2,1,1)
hist(modLow,[-1:0.1:1])
title('Modulation For Low Power, 500ms After Laser Offset')
subplot(2,1,2)
hist(modHi,[-1:0.1:1])
title('Modulation For High Power, 500ms After Laser Offset')

hFig = figure;
subplot(2,1,1)
hold on
for i = 1:length(msns)
    plot([testPre(i,1),testPost(i,1)],'c')
end
plot([mean(testPre(:,1)),mean(testPost(:,1))],'k','LineWidth',2)
xlim([0 3])

subplot(2,1,2)
hold on
for i = 1:length(msns)
    plot([testPre(i,2),testPost(i,2)],'c')
end
plot([mean(testPre(:,2)),mean(testPost(:,2))],'k','LineWidth',2)
xlim([0 3])


%okay, now lets try and bring out data regarding changes in firing rates
%over time. 


masterInd = 1;
lowInd = 1;
hiInd = 1;
for i = 1:numFiles
   %load the target file
    load(targetFiles{i})
    numUnits = length(s.DesignationName);
    %find time of first/last stimulations of low/hi power
    timePoints = s.LaserData.LaserStartTimes([1,100,101,200]);
    
    for j = 1:numUnits
        lowSpikes = s.(s.DesignationName{j}).SpikeTimes(s.(s.DesignationName{j}).SpikeTimes >  timePoints(1) & s.(s.DesignationName{j}).SpikeTimes< timePoints(2))-timePoints(1);
        hiSpikes = s.(s.DesignationName{j}).SpikeTimes(s.(s.DesignationName{j}).SpikeTimes >  timePoints(3) & s.(s.DesignationName{j}).SpikeTimes< timePoints(4))-timePoints(3);
        
        lowPowerStore(lowInd:lowInd + length(lowSpikes) -1) = lowSpikes;
        hiPowerStore(lowInd:lowInd + length(hiSpikes) -1) = hiSpikes;
        
        lowInd = lowInd + length(lowSpikes);
        hiInd = hiInd + length(hiSpikes);
    end
    masterInd = masterInd + numUnits;
    
    
    
end

%now try to do this while separating PV and MSN

masterInd = 1;
lowInd = 1;
hiInd = 1;
for i = 1:numFiles
   %load the target file
    load(targetFiles{i})
    msns =find(s.MasterSheet(:,1) ==0);
    pvs = find(s.MasterSheet(:,1) == 1);
    unitsMSN = s.DesignationName(msns);
    unitsPV = s.DesignationName(pvs);
    numUnits = length(unitsMSN);
%     numUnits = length(unitsPV);
%     numUnits = length(s.DesignationName);
    %find time of first/last stimulations of low/hi power
    timePoints = s.LaserData.LaserStartTimes([1,100,101,200]);
    
    for j = 1:numUnits
        lowSpikes = s.(unitsMSN{j}).SpikeTimes(s.(unitsMSN{j}).SpikeTimes >  timePoints(1) & s.(unitsMSN{j}).SpikeTimes< timePoints(2))-timePoints(1);
        hiSpikes = s.(unitsMSN{j}).SpikeTimes(s.(unitsMSN{j}).SpikeTimes >  timePoints(3) & s.(unitsMSN{j}).SpikeTimes< timePoints(4))-timePoints(3);
        
        lowPowerStore(lowInd:lowInd + length(lowSpikes) -1) = lowSpikes;
        hiPowerStore(lowInd:lowInd + length(hiSpikes) -1) = hiSpikes;
        
        lowInd = lowInd + length(lowSpikes);
        hiInd = hiInd + length(hiSpikes);
    end
    masterInd = masterInd + numUnits;
    
    
    
end




%now lets try and pull histograms in bins of 25 trials (4 per set).

%now lets pull out the average histogram
masterInd = 1;
for i = 1:numFiles
   %load the target file
    load(targetFiles{i})
    numUnits = length(s.DesignationName);
%     normHistStore = zeros(120,1,2);
%     zHistStore = zeros(120,1,2);
    for j = 1:numUnits
        normHistStore(:,masterInd + j - 1,1) = mean(s.(s.DesignationName{j}).TrialHists(:,1:25)');
        normHistStore(:,masterInd + j - 1,2) = mean(s.(s.DesignationName{j}).TrialHists(:,26:50)');
        normHistStore(:,masterInd + j - 1,3) = mean(s.(s.DesignationName{j}).TrialHists(:,51:75)');
        normHistStore(:,masterInd + j - 1,4) = mean(s.(s.DesignationName{j}).TrialHists(:,76:100)');
        normHistStore(:,masterInd + j - 1,5) = mean(s.(s.DesignationName{j}).TrialHists(:,101:125)');
        normHistStore(:,masterInd + j - 1,6) = mean(s.(s.DesignationName{j}).TrialHists(:,126:150)');
        normHistStore(:,masterInd + j - 1,7) = mean(s.(s.DesignationName{j}).TrialHists(:,151:175)');
        normHistStore(:,masterInd + j - 1,8) = mean(s.(s.DesignationName{j}).TrialHists(:,176:200)');
        infoStore(masterInd + j - 1) = s.MasterSheet(j,1);
    end
    masterInd = masterInd + numUnits;
    
    
    
end

testPVs = find(infoStore == 1);
testMSNs = find(infoStore == 0);

meanPVs = squeeze(mean(normHistStore(:,testPVs,:),2));
meanMSNs = squeeze(mean(normHistStore(:,testMSNs,:),2));

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure

subplot(2,1,1)
hold on
for i = 1:8
plot([1+(120*(i-1)):i*120],meanPVs(:,i)*20)
plot([i*120 i*120],[0 8],'k')
plot([1+(120*(i-1)),i*120],[mean(meanPVs(:,i)*20,1) mean(meanPVs(:,i)*20,1)],'r')
end

subplot(2,1,2)
hold on
for i = 1:8
plot([1+(120*(i-1)):i*120],meanMSNs(:,i)*20)
plot([i*120 i*120],[0 3],'k')
plot([1+(120*(i-1)),i*120],[mean(meanMSNs(:,i)*20,1) mean(meanMSNs(:,i)*20,1)],'r')
end



%try normalizing by baseline rates

%now lets pull out the average histogram
masterInd = 1;
for i = 1:numFiles
   %load the target file
    load(targetFiles{i})
    numUnits = length(s.DesignationName);
%     normHistStore = zeros(120,1,2);
%     zHistStore = zeros(120,1,2);
    for j = 1:numUnits
        baseRate = mean(s.(s.DesignationName{j}).HistogramLaser(1:20));
        normHistStore(:,masterInd + j - 1,1) = mean(s.(s.DesignationName{j}).TrialHists(:,1:25)')/baseRate;
        normHistStore(:,masterInd + j - 1,2) = mean(s.(s.DesignationName{j}).TrialHists(:,26:50)')/baseRate;
        normHistStore(:,masterInd + j - 1,3) = mean(s.(s.DesignationName{j}).TrialHists(:,51:75)')/baseRate;
        normHistStore(:,masterInd + j - 1,4) = mean(s.(s.DesignationName{j}).TrialHists(:,76:100)')/baseRate;
        normHistStore(:,masterInd + j - 1,5) = mean(s.(s.DesignationName{j}).TrialHists(:,101:125)')/baseRate;
        normHistStore(:,masterInd + j - 1,6) = mean(s.(s.DesignationName{j}).TrialHists(:,126:150)')/baseRate;
        normHistStore(:,masterInd + j - 1,7) = mean(s.(s.DesignationName{j}).TrialHists(:,151:175)')/baseRate;
        normHistStore(:,masterInd + j - 1,8) = mean(s.(s.DesignationName{j}).TrialHists(:,176:200)')/baseRate;
    end
    masterInd = masterInd + numUnits;
    
    
    
end

meanPVs = squeeze(mean(normHistStore(:,pvs,:),2));
meanMSNs = squeeze(mean(normHistStore(:,msns,:),2));


subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure

subplot(2,1,1)
hold on
for i = 1:8
plot([1+(120*(i-1)):i*120],meanPVs(:,i)*20)
plot([i*120 i*120],[0 2],'k')
plot([1+(120*(i-1)),i*120],[mean(meanPVs(:,i)*20,1) mean(meanPVs(:,i)*20,1)],'r')
end

subplot(2,1,2)
hold on
for i = 1:8
plot([1+(120*(i-1)):i*120],meanMSNs(:,i)*20)
plot([i*120 i*120],[0 3],'k')
plot([1+(120*(i-1)),i*120],[mean(meanMSNs(:,i)*20,1) mean(meanMSNs(:,i)*20,1)],'r')
end
