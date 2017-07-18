% testerWrapperWhiteLaser

%pull the analysis files that have already been processed.
% targetFiles = {'170629_ML170621G_R10_2050_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat',...
%     '170629_ML170621G_R10_2300_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat',...
%     '170629_ML170621G_R10_2550_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat',...
%     '170629_ML170621G_R10_2800_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat'};
% 
% targetFiles = {'170705_ML170621I_L10_2000_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat',...
%     '170705_ML170621I_L10_2250_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat',...
%     '170705_ML170621I_L10_2500_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat'};

targetFiles = {'170706_ML170621I_R10_2000_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat',...
    '170706_ML170621I_R10_2250_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat',...
    '170706_ML170621I_R10_2500_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis.mat'};

% 170705_ML170621I_L10_2000_3mWLaserWhiteNoiseThreepeatWhiteLaserComboAnalysis


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
    [indLaserOnly] = functionCellStringFind(masterHeader,'LaserOnlyAverage');
    [indToneOnly] = functionCellStringFind(masterHeader,'ToneOnlyAverage');
    [indToneLaser] = functionCellStringFind(masterHeader,'ToneLaserAverage');
    [indPreAverageRunning] = functionCellStringFind(masterHeader,'Pre AverageRunning');
    [indPreAverageStationary] = functionCellStringFind(masterHeader,'Pre AverageStationary');
    [indLocoAUC] = functionCellStringFind(masterHeader,'LocoAUCScore');
    [indLocoSig] = functionCellStringFind(masterHeader,'LocoAUCSignificance');
    
    %now i need to adjust the depth system. First, add the height of the
    %probe in terms of number of probe sites. then multiply by negative 1 and 12
    master(:,indDistance) = (master(:,indDistance) + s.ShankLength)*-12.5 + probeDepth(i);
    numUnits = size(master,1);
    fullMaster(masterInd: masterInd + numUnits - 1,:) = master;
    fullMasterHeaders(:,i) = masterHeader;
    masterInd = masterInd + numUnits;
    
    
end


%% This is some basic analysis for looking at things 170705
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%calculate modulation ind for laser only and tone only

mod1 = (fullMaster(:,19)-fullMaster(:,15))./(fullMaster(:,19)+fullMaster(:,15));

%alternative mod1, which also factors in post?
tester = (fullMaster(:,19)-fullMaster(:,14))./(fullMaster(:,19)+fullMaster(:,14));

mod2 = (fullMaster(:,21)-fullMaster(:,20))./(fullMaster(:,21)+fullMaster(:,20));
scottMod1 = (fullMaster(:,19))./(fullMaster(:,19)+fullMaster(:,15));
scottMod2 = (fullMaster(:,21))./(fullMaster(:,21)+fullMaster(:,20));

hFig = figure

%find pvs and msn
msns = find(fullMaster(:,1) ==0);
pvs = find(fullMaster(:,1) ==1);

%column 1
subplot(4,3,1)
plot(fullMaster(:,2),fullMaster(:,3),'b.')
hold on
plot(fullMaster(pvs,2),fullMaster(pvs,3),'ro')
title('PeakTrough(x) vs Spike Rate(y), PV in red')

subplot(4,3,4)
hist(fullMaster(pvs,3),100)
title('Histogram of PV Firing Rate')

subplot(4,3,7)
hist(fullMaster(msns,3),100)
title('Histogram of MSN Firing Rate')

subplot(4,3,10)
plot(fullMaster(msns,3),mod1(msns),'b.')
hold on
plot([0 max(fullMaster(msns,3))],[0 0],'k')
ylim([-1 1])
title('MSN Firing Rate vs Modulation Index')

%column two
subplot(4,6,3)
%find pv neurons on shank 1
pvShank1 = find(fullMaster(:,1) == 1 & fullMaster(:,5) == 1);
pvShank2 = find(fullMaster(:,1) == 1 & fullMaster(:,5) == 2);

for i = 1:length(pvShank1)
    hold on
    plot([0 mod1(pvShank1(i))],[fullMaster(pvShank1(i),4) fullMaster(pvShank1(i),4)],'b')
    plot(mod1(pvShank1(i)),fullMaster(pvShank1(i),4),'b.')
    plot([0 mod2(pvShank1(i))],[fullMaster(pvShank1(i),4) fullMaster(pvShank1(i),4)],'g')
    plot(mod2(pvShank1(i)),fullMaster(pvShank1(i),4),'g.')
end

plot([0 0],[min(fullMaster(:,4)),max(fullMaster(:,4))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,4))-10,max(fullMaster(:,4))+10])

title('PV Mod,Laser(b) Tone(g) by Depth Shank 1')

subplot(4,6,4)

for i = 1:length(pvShank2)
    hold on
    plot([0 mod1(pvShank2(i))],[fullMaster(pvShank2(i),4) fullMaster(pvShank2(i),4)],'b')
    plot(mod1(pvShank2(i)),fullMaster(pvShank2(i),4),'b.')
    plot([0 mod2(pvShank2(i))],[fullMaster(pvShank2(i),4) fullMaster(pvShank2(i),4)],'g')
    plot(mod2(pvShank2(i)),fullMaster(pvShank2(i),4),'g.')
end

plot([0 0],[min(fullMaster(:,4)),max(fullMaster(:,4))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,4))-10,max(fullMaster(:,4))+10])

title('PV Mod,Laser(b) Tone(g) by Depth Shank 2')


subplot(4,6,9)
hist(mod1(pvs),[-1:0.1:1])
xlim([-1 1])
title('MODLaser')


subplot(4,6,10)
hist(mod2(pvs),[-1:0.1:1])
xlim([-1 1])
title('MODTone')



%column three
subplot(4,6,5)
%find msn neurons on shank 1
msnShank1 = find(fullMaster(:,1) == 0 & fullMaster(:,5) == 1);
msnShank2 = find(fullMaster(:,1) == 0 & fullMaster(:,5) == 2);

for i = 1:length(msnShank1)
    hold on
    plot([0 mod1(msnShank1(i))],[fullMaster(msnShank1(i),4) fullMaster(msnShank1(i),4)],'b')
    plot(mod1(msnShank1(i)),fullMaster(msnShank1(i),4),'b.')
    plot([0 mod2(msnShank1(i))],[fullMaster(msnShank1(i),4) fullMaster(msnShank1(i),4)],'g')
    plot(mod2(msnShank1(i)),fullMaster(msnShank1(i),4),'g.')
end

plot([0 0],[min(fullMaster(:,4)),max(fullMaster(:,4))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,4))-10,max(fullMaster(:,4))+10])

title('MSN Mod,Laser(b) Tone(g) by Depth Shank 1')

subplot(4,6,6)

for i = 1:length(msnShank2)
    hold on
    plot([0 mod1(msnShank2(i))],[fullMaster(msnShank2(i),4) fullMaster(msnShank2(i),4)],'b')
    plot(mod1(msnShank2(i)),fullMaster(msnShank2(i),4),'b.')
    plot([0 mod2(msnShank2(i))],[fullMaster(msnShank2(i),4) fullMaster(msnShank2(i),4)],'g')
    plot(mod2(msnShank2(i)),fullMaster(msnShank2(i),4),'g.')
end

plot([0 0],[min(fullMaster(:,4)),max(fullMaster(:,4))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,4))-10,max(fullMaster(:,4))+10])

title('MSN Mod,Laser(b) Tone(g) by Depth Shank 2')


subplot(4,6,11)
hist(mod1(msns),[-1:0.1:1])
xlim([-1 1])
title('MODLaser')


subplot(4,6,12)
hist(mod2(msns),[-1:0.1:1])
xlim([-1 1])
title('MODTone')

%find only tone responsive msns
msnToneShank1 = find(fullMaster(:,1) == 0 & fullMaster(:,5) == 1 & fullMaster(:,6) == 1);
msnToneShank2 = find(fullMaster(:,1) == 0 & fullMaster(:,5) == 2 & fullMaster(:,6) == 1);
msnTones = find(fullMaster(:,1) == 0 & fullMaster(:,6) == 1);

subplot(4,6,17)

for i = 1:length(msnToneShank1)
    hold on
    plot([0 mod1(msnToneShank1(i))],[fullMaster(msnToneShank1(i),4) fullMaster(msnToneShank1(i),4)],'b')
    plot(mod1(msnToneShank1(i)),fullMaster(msnToneShank1(i),4),'b.')
    plot([0 mod2(msnToneShank1(i))],[fullMaster(msnToneShank1(i),4) fullMaster(msnToneShank1(i),4)],'g')
    plot(mod2(msnToneShank1(i)),fullMaster(msnToneShank1(i),4),'g.')
end

plot([0 0],[min(fullMaster(:,4)),max(fullMaster(:,4))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,4))-10,max(fullMaster(:,4))+10])

title('MSNTONE Mod,Laser(b) Tone(g) by Depth Shank 1')

subplot(4,6,18)

for i = 1:length(msnToneShank2)
    hold on
    plot([0 mod1(msnToneShank2(i))],[fullMaster(msnToneShank2(i),4) fullMaster(msnToneShank2(i),4)],'b')
    plot(mod1(msnToneShank2(i)),fullMaster(msnToneShank2(i),4),'b.')
    plot([0 mod2(msnToneShank2(i))],[fullMaster(msnToneShank2(i),4) fullMaster(msnToneShank2(i),4)],'g')
    plot(mod2(msnToneShank2(i)),fullMaster(msnToneShank2(i),4),'g.')
end

plot([0 0],[min(fullMaster(:,4)),max(fullMaster(:,4))],'k')
set (gca,'Ydir','reverse')
xlim([-1 1])
ylim([min(fullMaster(:,4))-10,max(fullMaster(:,4))+10])

title('MSNTONE Mod,Laser(b) Tone(g) by Depth Shank 2')


subplot(4,6,23)
hist(mod1(msnTones),[-1:0.1:1])
xlim([-1 1])
title('MODLaser')


subplot(4,6,24)
hist(mod2(msnTones),[-1:0.1:1])
xlim([-1 1])
title('MODTone')


