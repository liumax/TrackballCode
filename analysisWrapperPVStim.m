% testerWrapperWhiteLaser

% % % %pull the analysis files that have already been processed.
% targetFiles = {'170714_ML170620C_L10_2000_3mW2sOn4sOffLaserStimAnalysis.mat',...
%     '170714_ML170620C_L10_2250_3mW2sOn4sOffLaserStimAnalysis.mat',...
%     '170714_ML170620C_L10_2500_3mW2sOn4sOffLaserStimAnalysis.mat'};

% targetFiles = {'170714_ML170620C_R10_2000_3mW2sOn4sOffLaserStimAnalysis.mat',...
%     '170714_ML170620C_R10_2250_Second3mW2sOn4sOffLaserStimAnalysis.mat'};

%get the number of files





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
    masterHeader{end+1} = 'SpikeWidth';
    masterHeader{end+1} = 'signRankSigOfRestricted';
    masterHeader{end+1} = 'signRankSigOfNorm';
    disp(strcat('Headings Generated for File',num2str(i),'/',num2str(numFiles)))
    %170905 need to do some repairs to add percentage of trials with
    %firing. 
    numCells = length(s.DesignationName);
    numStims = length(s.Timing.LaserTimes);
    laserRasterData = zeros(numCells,1);
    isiRatio = zeros(numCells,1);
    spikeWidth = zeros(numCells,1);
    resSig = zeros(numCells,1);
    normSig = zeros(numCells,1);
    numSpikes = zeros(numCells,1);
    preITIratio = zeros(numCells,1);
    preITInum = zeros(numCells,1);
    isiCov = zeros(numCells,1);
    for j = 1:numCells
        disp(strcat('Analyzing Cell #',num2str(j),'/',num2str(numCells)))
        %pull rasterLaser
        laserRasterData(j) = length(unique(s.(s.DesignationName{j}).RasterLaser(:,2)));
        %also pull isi data. This is crude measure over entire firing
        %period. 
        isiTimes = diff(s.(s.DesignationName{j}).SpikeTimes);
        isiCov(j) = std(isiTimes)/mean(isiTimes);
        isiRatio(j) = length(find(isiTimes < 0.033 & isiTimes > 0.01))/length(isiTimes);
        numSpikes(j) = length(s.(s.DesignationName{j}).SpikeTimes);
        %now lets try pulling just from the baseline
        itiCounter = 1;
        preITIstore = [];
        for k = 1:length(s.LaserData.LaserStartTimes)
            %first, pull out spikes from ITIs
            if k == 1
                selSpikes = s.(s.DesignationName{j}).SpikeTimes(s.(s.DesignationName{j}).SpikeTimes < s.LaserData.LaserStartTimes(k));
            else
                selSpikes = s.(s.DesignationName{j}).SpikeTimes(s.(s.DesignationName{j}).SpikeTimes < s.LaserData.LaserStartTimes(k));
                selSpikes = selSpikes(selSpikes>(s.LaserData.LaserEndTimes(k-1)+1));%adding in 1 second fudge factor to avoid the period immediately following laser offset. 
            end
            lengthSpikes = length(selSpikes);
            if lengthSpikes > 1
                preDiff = diff(selSpikes);
                preITIstore(itiCounter:itiCounter + length(preDiff)-1) = preDiff;
                itiCounter = itiCounter + length(preDiff);
            end
        end
        preITIratio(j) = length(find(preITIstore < 0.033 & preITIstore > 0.01))/length(preITIstore);
        preITInum(j) = length(preITIstore);
        %now lets try pulling spike width (half max)
        waves = s.(s.DesignationName{j}).AverageWaveForms;
        maxWave = max(waves);
        [maxVal maxInd] = max(maxWave);
        %chose the big wave, interpolate to fine degree
        chosenWave = waves(:,maxInd);
        interpVect = [1:0.1:40];
        interpWave = interp1(1:40,chosenWave,interpVect,'spline');
        
        [pkVal pkInd] = max(interpWave(100:end));
        pkInd = pkInd + 100 - 1;
        %now we need to crawl forward and backwards to find the first point
        %equal to half  the peak value
        firstHalf = find(interpWave(1:pkInd) <= pkVal/2,1,'last');
        secondHalf = find(interpWave(pkInd:end) >= pkVal/2,1,'last') + pkInd;
        spikeWidth(j) = (secondHalf - firstHalf)/300000;
        %now lets try calculating p value using signed rank.
        avFire = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(:,[4,5])');
        laserFire = s.(s.DesignationName{j}).TrialBinnedSpikes(:,6);
        resSig(j) = signrank(avFire,laserFire);
        avFire = mean(s.(s.DesignationName{j}).TrialBinnedSpikes(:,[1,2])');
        laserFire = s.(s.DesignationName{j}).TrialBinnedSpikes(:,3);
        normSig(j) = signrank(avFire,laserFire);
        
    end
    disp('Finished Cell Data, Storing in Master')
    laserRasterData = laserRasterData/numStims;
    
    
    
    %now i need to adjust the depth system. First, add the height of the
    %probe in terms of number of probe sites. then multiply by negative 1 and 12
    master(:,indDistance) = (master(:,indDistance) + s.ShankLength)*-12.5 + probeDepth(i);
    numUnits = size(master,1);
    fullMaster(masterInd: masterInd + numUnits - 1,1:size(master,2)) = master;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+1) = laserRasterData;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+2) = numStims;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+3) = isiRatio;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+4) = spikeWidth;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+5) = resSig;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+6) = normSig;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+7) = numSpikes;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+8) = preITIratio;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+9) = preITInum;
    fullMaster(masterInd: masterInd + numUnits - 1,size(master,2)+10) = isiCov;
    fullMasterHeaders(:,i) = masterHeader;
    masterInd = masterInd + numUnits;
    
    
end


%% This is some basic analysis for looking at things 170705
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%calculate modulation ind for normal and restricted periods
mod1 = (fullMaster(:,indLaserAverage)-fullMaster(:,indPreAverage))./(fullMaster(:,indLaserAverage)+fullMaster(:,indPreAverage));
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
hist(fullMaster(:,20),[0:0.01:1])
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

