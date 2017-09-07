% testerWrapperWhiteLaser

% % % %pull the analysis files that have already been processed.
% targetFiles = {'170714_ML170620C_L10_2000_3mW2sOn4sOffLaserStimAnalysis.mat',...
%     '170714_ML170620C_L10_2250_3mW2sOn4sOffLaserStimAnalysis.mat',...
%     '170714_ML170620C_L10_2500_3mW2sOn4sOffLaserStimAnalysis.mat'};

% targetFiles = {'170714_ML170620C_R10_2000_3mW2sOn4sOffLaserStimAnalysis.mat',...
%     '170714_ML170620C_R10_2250_Second3mW2sOn4sOffLaserStimAnalysis.mat'};

%get the number of files
numFiles = length(targetFiles);

%extract distances
for i = 1:numFiles
    targetName = targetFiles{i};
    findSpace = strfind(targetName,'_');
    probeDepth(i) = str2num(targetName(findSpace(3)+1:findSpace(4)-1));
end

%now we need to make a big master sheet.
masterInd = 1;
for i = 1:numFiles
    %load the target file
    load(targetFiles{i})
    %extract master array
    master = s.MasterSheet;
    masterHeader = s.MasterDesigs;
    %pull indices
    [indPVMSN] = functionCellStringFind(masterHeader,'PV/MSN Desig');
    [indPkTrough] = functionCellStringFind(masterHeader,'PeakTroughTimeinMS');
    [indOverFire] = functionCellStringFind(masterHeader,'OverallFiringRate');
    [indDistance] = functionCellStringFind(masterHeader,'Distance from Top of Shank');
    [indShankDes] = functionCellStringFind(masterHeader,'Shank Designation');
    [indPreAverage] = functionCellStringFind(masterHeader,'PreAverage');
    [indPostAverage] = functionCellStringFind(masterHeader,'PostAverage');
    [indLaserAverage] = functionCellStringFind(masterHeader,'LaserAverage');
    [indPreResAverage] = functionCellStringFind(masterHeader,'PreResAverage');
    [indPostResAverage] = functionCellStringFind(masterHeader,'PostResAverage');
    [indLaserResAverage] = functionCellStringFind(masterHeader,'LaserResAverage');
    [indLocoAUC] = functionCellStringFind(masterHeader,'LocoAUCScore');
    [indLocoSig] = functionCellStringFind(masterHeader,'LocoAUCSignificance');
    masterHeader{end+1} = 'PercentFiring';
    masterHeader{end+1} = 'NumberStimulations';
    %170905 need to do some repairs to add percentage of trials with
    %firing. 
    numCells = length(s.DesignationName);
    numStims = length(s.Timing.LaserTimes);
    laserRasterData = zeros(numCells,1);
    for j = 1:numCells
        %pull rasterLaser
        laserRasterData(j) = length(unique(s.(s.DesignationName{j}).RasterLaser(:,2)));
    end
    
    laserRasterData = laserRasterData/numStims;
    
    
    
    %now i need to adjust the depth system. First, add the height of the
    %probe in terms of number of probe sites. then multiply by negative 1 and 12
    master(:,indDistance) = (master(:,indDistance) + s.ShankLength)*-12.5 + probeDepth(i);
    numUnits = size(master,1);
    fullMaster(masterInd: masterInd + numUnits - 1,1:size(master,2)) = master;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+1) = laserRasterData;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+2) = numStims;
    fullMasterHeaders(:,i) = masterHeader;
    masterInd = masterInd + numUnits;
    
    
end


%% This is some basic analysis for looking at things 170705
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%calculate modulation ind for normal and restricted periods
mod1 = (fullMaster(:,indLaserAverage)-fullMaster(:,indPreAverage))./(fullMaster(:,indLaserResAverage)+fullMaster(:,indPreAverage));
scottMod1 = (fullMaster(:,indLaserAverage))./(fullMaster(:,indLaserAverage)+fullMaster(:,indPreAverage));
%find for restricted
modRes1 = (fullMaster(:,indLaserResAverage)-fullMaster(:,indPreResAverage))./(fullMaster(:,indLaserResAverage)+fullMaster(:,indPreResAverage));
scottModRes1 = (fullMaster(:,indLaserResAverage))./(fullMaster(:,indLaserResAverage)+fullMaster(:,indPreResAverage));

hFig = figure
set(hFig, 'Position', [10 80 1240 850])

%find pvs and msn
msns = find(fullMaster(:,indPVMSN) ==0);
pvs = find(fullMaster(:,indPVMSN) ==1);

%column 1
subplot(4,3,1)
plot(fullMaster(:,indPkTrough),fullMaster(:,indOverFire),'b.')
hold on
plot(fullMaster(pvs,indPkTrough),fullMaster(pvs,indOverFire),'ro')
title('PeakTrough(x) vs Spike Rate(y), PV in red')

subplot(4,3,4)
hist(fullMaster(:,end-1),[0:0.01:1])
title('Histogram of % Spiking Per Trial')



%column two
subplot(4,6,3)
%find pv neurons on shank 1
pvShank1 = find(fullMaster(:,indPVMSN) == 1 & fullMaster(:,indShankDes) == 1);
pvShank2 = find(fullMaster(:,indPVMSN) == 1 & fullMaster(:,indShankDes) == 2);

for i = 1:length(pvShank1)
    hold on
    plot([0 mod1(pvShank1(i))],[fullMaster(pvShank1(i),indDistance) fullMaster(pvShank1(i),indDistance)],'b')
    plot(mod1(pvShank1(i)),fullMaster(pvShank1(i),indDistance),'b.')
end

plot([0 0],[min(fullMaster(:,indDistance)),max(fullMaster(:,indDistance))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,indDistance))-10,max(fullMaster(:,indDistance))+10])

title('PV Mod Shank 1')

subplot(4,6,4)

for i = 1:length(pvShank2)
    hold on
    plot([0 mod1(pvShank2(i))],[fullMaster(pvShank2(i),indDistance) fullMaster(pvShank2(i),indDistance)],'b')
    plot(mod1(pvShank2(i)),fullMaster(pvShank2(i),indDistance),'b.')
end

plot([0 0],[min(fullMaster(:,indDistance)),max(fullMaster(:,indDistance))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,indDistance))-10,max(fullMaster(:,indDistance))+10])

title('PV Mod Shank 2')


subplot(4,6,9)
hist(mod1(pvs),[-1:0.1:1])
xlim([-1 1])
title('MODLaser')


subplot(4,6,10)
hist(modRes1(pvs),[-1:0.1:1])
xlim([-1 1])
title('MODRestricted')


subplot(4,3,8)
hist(fullMaster(pvs,indOverFire),100)
title(strcat('Histogram of PV Firing Rate',num2str(length(pvs))))

subplot(4,3,11)
plot(fullMaster(pvs,indOverFire),mod1(pvs),'b.')
hold on
plot([0 max(fullMaster(pvs,indOverFire))],[0 0],'k')
ylim([-1 1])
title('PV Firing Rate vs Modulation Index')

%column three
subplot(4,6,5)
%find msn neurons on shank 1
msnShank1 = find(fullMaster(:,indPVMSN) == 0 & fullMaster(:,indShankDes) == 1);
msnShank2 = find(fullMaster(:,indPVMSN) == 0 & fullMaster(:,indShankDes) == 2);

for i = 1:length(msnShank1)
    hold on
    plot([0 mod1(msnShank1(i))],[fullMaster(msnShank1(i),indDistance) fullMaster(msnShank1(i),indDistance)],'b')
    plot(mod1(msnShank1(i)),fullMaster(msnShank1(i),indDistance),'b.')
end

plot([0 0],[min(fullMaster(:,indDistance)),max(fullMaster(:,indDistance))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,indDistance))-10,max(fullMaster(:,indDistance))+10])

title('MSN Mod Shank 1')

subplot(4,6,6)

for i = 1:length(msnShank2)
    hold on
    plot([0 mod1(msnShank2(i))],[fullMaster(msnShank2(i),indDistance) fullMaster(msnShank2(i),indDistance)],'b')
    plot(mod1(msnShank2(i)),fullMaster(msnShank2(i),indDistance),'b.')
end

plot([0 0],[min(fullMaster(:,indDistance)),max(fullMaster(:,indDistance))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,indDistance))-10,max(fullMaster(:,indDistance))+10])

title('MSN Mod Shank 2')


subplot(4,6,11)
hist(mod1(msns),[-1:0.1:1])
xlim([-1 1])
title('MODLaser')


subplot(4,6,12)
hist(modRes1(msns),[-1:0.1:1])
xlim([-1 1])
title('MODRestricted')


subplot(4,3,9)
hist(fullMaster(msns,indOverFire),100)
title(strcat('Histogram of MSN Firing Rate',num2str(length(msns))))
% title('Histogram of MSN Firing Rate')

subplot(4,3,12)
plot(fullMaster(msns,indOverFire),mod1(msns),'b.')
hold on
plot([0 max(fullMaster(msns,indOverFire))],[0 0],'k')
ylim([-1 1])
title('MSN Firing Rate vs Modulation Index')


%plot ROC vs modulation

subplot(2,3,4)
hold on
plot(mod1(msns),fullMaster(msns,indLocoAUC),'k.')
plot([0 0],[0 1])
plot([-1 1],[0.5 0.5])
findSig = find(fullMaster(:,indLocoSig) == 1);
findSig = intersect(msns,findSig);
plot(mod1(findSig),fullMaster(findSig,indLocoAUC),'ro')
xlim([-1 1])
ylim([0 1])
title('Mod Index (x) vs Loco AUC (R = sig)')


spikeGraphName = 'LaserStimSummary';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


save('AnalysisResults','fullMaster','fullMasterHeaders')

