%this is test code to pull out 50 pulse laser high low data

fileNames = what;
fileNames = fileNames.mat;
counter = 1;
for i = 1:length(fileNames)
    load(fileNames{i})
    
    %pull master sheet
    masterDesig = s.MasterDesigs;
    masterSheet(counter:counter+size(s.MasterSheet,1)-1,:) = s.MasterSheet;
    countSheet(counter:counter+size(s.MasterSheet,1)-1) = i;
    %pull average rates
    for j = 1:length(s.DesignationName)
        averageHistLow(counter-1+j,:) = s.(s.DesignationName{j}).LaserHistLow;
        averageHistHi(counter-1+j,:) = s.(s.DesignationName{j}).LaserHistHi;
        
        %now i need to pull just laser baseline shit
        lowTrials = [(s.LowBlocks(1)-1)*50+1:(s.LowBlocks(1)-1)*50+50,(s.LowBlocks(2)-1)*50+1:(s.LowBlocks(2)-1)*50+50];
        hiTrials = [(s.HiBlocks(1)-1)*50+1:(s.HiBlocks(1)-1)*50+50,(s.HiBlocks(2)-1)*50+1:(s.HiBlocks(2)-1)*50+50];
        rastersLow = [];
        rastersHi = [];
        rastersLow = s.(s.DesignationName{j}).RasterLaser(ismember(s.(s.DesignationName{j}).RasterLaser(:,2),lowTrials),1);
        rastersHi = s.(s.DesignationName{j}).RasterLaser(ismember(s.(s.DesignationName{j}).RasterLaser(:,2),hiTrials),1);
        %remove values greater than zero
        rastersLow(rastersLow>0) = [];
        rastersHi(rastersHi>0) = [];
        baselineVector = [-5:0.05:0];
        baselineLow = hist(rastersLow,baselineVector)/100/0.05;
        baselineHi = hist(rastersHi,baselineVector)/100/0.05;
        stdLow = std(baselineLow);
        stdHi = std(baselineHi);
        meanLow = mean(baselineLow);
        meanHi = mean(baselineHi);
        averageZLow(counter-1+j,:) = (s.(s.DesignationName{j}).LaserHistLow - meanLow)/stdLow;
        averageZHi(counter-1+j,:) = (s.(s.DesignationName{j}).LaserHistHi - meanHi)/stdHi;
        isiCov(counter - 1 + j) = std(diff(s.(s.DesignationName{j}).SpikeTimes))/mean(diff(s.(s.DesignationName{j}).SpikeTimes));
    end
        
    
    jumpsBack = round((s.Parameters.RasterWindow(1))/s.Parameters.InterpolationStepRotary);
    jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);
    
    
    velRaster = s.VelocityRasterLow;
    %pull velocity rasters
    
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
    fullVelFind = [1:length(s.Timing.LaserTimes)];
    fullVelFind(combFind) = [];
    
    typeStoreLow(i,:)  = [length(fullVelFind) length(fullStopZeroFind) length(laserZeroFind) length(preZeroFind) length(s.Timing.LaserTimes)];
    
    
    velRaster = s.VelocityRasterHi;
    %pull velocity rasters
    
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
    fullVelFind = [1:length(s.Timing.LaserTimes)];
    fullVelFind(combFind) = [];
    
    typeStoreHi(i,:)  = [length(fullVelFind) length(fullStopZeroFind) length(laserZeroFind) length(preZeroFind) length(s.Timing.LaserTimes)];
    
    counter = counter + size(s.MasterSheet,1);
end


pvs = find(masterSheet(:,1) == 1 & isiCov' > 1.1);
msns = find(masterSheet(:,1) == 0 & isiCov' > 1.1);

msnLowRateMean = mean(averageHistLow(msns,:));
msnHiRateMean = mean(averageHistHi(msns,:));

msnLowZMean = mean(averageZLow(msns,:));
msnHiZMean = mean(averageZHi(msns,:));

pvLowRateMean = mean(averageHistLow(pvs,:));
pvHiRateMean = mean(averageHistHi(pvs,:));

pvLowZMean = mean(averageZLow(pvs,:));
pvHiZMean = mean(averageZHi(pvs,:));

% msnLowRateMean = interp1([1:300],msnLowRateMean,[1:5/2:300]);
% msnHiRateMean = interp1([1:300],msnHiRateMean,[1:5/2:300]);
% 
% msnLowZMean = interp1([1:300],msnLowZMean,[1:5/2:300]);
% msnHiZMean = interp1([1:300],msnHiZMean,[1:5/2:300]);
% 
% pvLowRateMean = interp1([1:300],pvLowRateMean,[1:5/2:300]);
% pvHiRateMean = interp1([1:300],pvHiRateMean,[1:5/2:300]);
% 
% pvLowZMean = interp1([1:300],pvLowZMean,[1:5/2:300]);
% pvHiZMean = interp1([1:300],pvHiZMean,[1:5/2:300]);


msnLowRateMean = interp1([1:300],msnLowRateMean,[1:4:300]);
msnHiRateMean = interp1([1:300],msnHiRateMean,[1:4:300]);

msnLowZMean = interp1([1:300],msnLowZMean,[1:4:300]);
msnHiZMean = interp1([1:300],msnHiZMean,[1:4:300]);

pvLowRateMean = interp1([1:300],pvLowRateMean,[1:4:300]);
pvHiRateMean = interp1([1:300],pvHiRateMean,[1:4:300]);

pvLowZMean = interp1([1:300],pvLowZMean,[1:4:300]);
pvHiZMean = interp1([1:300],pvHiZMean,[1:4:300]);

msnLowRateStd = std(averageHistLow(msns,:))/sqrt(length(msns));
msnHiRateStd = std(averageHistHi(msns,:))/sqrt(length(msns));

msnLowZStd = std(averageZLow(msns,:))/sqrt(length(msns));
msnHiZStd = std(averageZHi(msns,:))/sqrt(length(msns));

pvLowRateStd = std(averageHistLow(pvs,:))/sqrt(length(pvs));
pvHiRateStd = std(averageHistHi(pvs,:))/sqrt(length(pvs));

pvLowZStd = std(averageZLow(pvs,:))/sqrt(length(pvs));
pvHiZStd = std(averageZHi(pvs,:))/sqrt(length(pvs));

% msnLowRateStd = interp1([1:300],msnLowRateStd,[1:5/2:300]);
% msnHiRateStd = interp1([1:300],msnHiRateStd,[1:5/2:300]);
% 
% msnLowZStd = interp1([1:300],msnLowZStd,[1:5/2:300]);
% msnHiZStd = interp1([1:300],msnHiZStd,[1:5/2:300]);
% 
% pvLowRateStd = interp1([1:300],pvLowRateStd,[1:5/2:300]);
% pvHiRateStd = interp1([1:300],pvHiRateStd,[1:5/2:300]);
% 
% pvLowZStd = interp1([1:300],pvLowZStd,[1:5/2:300]);
% pvHiZStd = interp1([1:300],pvHiZStd,[1:5/2:300]);

msnLowRateStd = interp1([1:300],msnLowRateStd,[1:4:300]);
msnHiRateStd = interp1([1:300],msnHiRateStd,[1:4:300]);

msnLowZStd = interp1([1:300],msnLowZStd,[1:4:300]);
msnHiZStd = interp1([1:300],msnHiZStd,[1:4:300]);

pvLowRateStd = interp1([1:300],pvLowRateStd,[1:4:300]);
pvHiRateStd = interp1([1:300],pvHiRateStd,[1:4:300]);

pvLowZStd = interp1([1:300],pvLowZStd,[1:4:300]);
pvHiZStd = interp1([1:300],pvHiZStd,[1:4:300]);

%readjust binning so that matches number of points with other recordings.
%this should make for 120 bins



hFig = figure;
hold on
plot(msnLowRateMean,'b','LineWidth',2)
plot(msnLowRateMean-msnLowRateStd,'b','LineWidth',2)
plot(msnLowRateMean+msnLowRateStd,'b','LineWidth',2)

plot(msnHiRateMean,'g','LineWidth',2)
plot(msnHiRateMean-msnHiRateStd,'g','LineWidth',2)
plot(msnHiRateMean+msnHiRateStd,'g','LineWidth',2)

xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
title('MSN Low (b) vs Hi (g)')
set(gca,'TickDir','out');
savefig(hFig,'msnHzLowVsHi');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'msnHzLowVsHi','-dpdf','-r0')

hFig = figure;
hold on
plot(msnLowZMean,'b','LineWidth',2)
plot(msnLowZMean-msnLowZStd,'b','LineWidth',2)
plot(msnLowZMean+msnLowZStd,'b','LineWidth',2)
plot(msnHiZMean,'g','LineWidth',2)
plot(msnHiZMean-msnHiZStd,'g','LineWidth',2)
plot(msnHiZMean+msnHiZStd,'g','LineWidth',2)
set(gca,'TickDir','out');
xlabel('Time (s)')
ylabel('Firing Rate (Z)')
title('MSN Low (b) vs Hi (g) Z')
savefig(hFig,'msnZLowVsHi');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'msnZLowVsHi','-dpdf','-r0')

hFig = figure;
hold on
plot(pvLowRateMean,'b','LineWidth',2)
plot(pvLowRateMean-pvLowRateStd,'b','LineWidth',2)
plot(pvLowRateMean+pvLowRateStd,'b','LineWidth',2)
plot(pvHiRateMean,'g','LineWidth',2)
plot(pvHiRateMean-pvHiRateStd,'g','LineWidth',2)
plot(pvHiRateMean+pvHiRateStd,'g','LineWidth',2)
set(gca,'TickDir','out');
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
title('PV Low (b) vs Hi (g)')
savefig(hFig,'pvHzLowVsHi');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'pvHzLowVsHi','-dpdf','-r0')

hFig = figure;
hold on
plot(pvLowZMean,'b','LineWidth',2)
plot(pvLowZMean-pvLowZStd,'b','LineWidth',2)
plot(pvLowZMean+pvLowZStd,'b','LineWidth',2)
plot(pvHiZMean,'g','LineWidth',2)
plot(pvHiZMean-pvHiZStd,'g','LineWidth',2)
plot(pvHiZMean+pvHiZStd,'g','LineWidth',2)
set(gca,'TickDir','out');
xlabel('Time (s)')
ylabel('Firing Rate (Z)')
title('PV Low (b) vs Hi (g) Z')
savefig(hFig,'pvZLowVsHi');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'pvZLowVsHi','-dpdf','-r0')


%calculate modulation index
modLow = (masterSheet(:,20) - (mean(masterSheet(:,18:19)'))')./(masterSheet(:,20) + (mean(masterSheet(:,18:19)'))');
modHi = (masterSheet(:,23) - (mean(masterSheet(:,21:22)'))')./(masterSheet(:,23) + (mean(masterSheet(:,21:22)'))');

hFig = figure
subplot(2,1,1)
hist(modLow(msns),[-1:0.1:1])
xlim([-1 1])
ylim([0 30])
title('Modulation Index MSNs, Low Power (3 mW)')
set(gca,'TickDir','out');
subplot(2,1,2)
hist(modHi(msns),[-1:0.1:1])
xlim([-1 1])
ylim([0 30])
title('Modulation Index MSNs, Hi Power (15 mW)')
set(gca,'TickDir','out');
savefig(hFig,'msnModHist');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'msnModHist','-dpdf','-r0')

hFig = figure
subplot(2,1,1)
hist(modLow(pvs),[-1:0.1:1])
xlim([-1 1])
ylim([0 10])
title('Modulation Index PV, Low Power (3 mW)')
set(gca,'TickDir','out');
subplot(2,1,2)
hist(modHi(pvs),[-1:0.1:1])
xlim([-1 1])
ylim([0 10])
title('Modulation Index PV, Hi Power (15 mW)')
set(gca,'TickDir','out');
savefig(hFig,'pvModHist');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'pvModHist','-dpdf','-r0')

pvalmsn = signrank(modLow(msns),modHi(msns));
pvalpv = signrank(modLow(pvs),modHi(pvs));


hFig= figure;
subplot(2,1,1)
hold on
for i = 1:length(msns)
    plot([modLow(msns(i)) modHi(msns(i))],'k')
end
plot(ones(length(msns),1),modLow(msns),'b.')
plot(ones(length(msns),1)*2,modHi(msns),'g.')
errorbar([mean(modLow(msns)) mean(modHi(msns))],[std(modLow(msns))/sqrt(length(msns)) std(modHi(msns))/sqrt(length(msns))],'LineWidth',2)
xlim([0.8 2.2])
ylim([-1 1])
xlabel('3 mW vs 15 mW')
ylabel('Modulation Index')
set(gca,'TickDir','out');
title(strcat('msns Mod Index. Pval:',num2str(pvalmsn)))

subplot(2,1,2)
hold on
for i = 1:length(pvs)
    plot([modLow(pvs(i)) modHi(pvs(i))],'k')
end
plot(ones(length(pvs),1),modLow(pvs),'b.')
plot(ones(length(pvs),1)*2,modHi(pvs),'g.')
errorbar([mean(modLow(pvs)) mean(modHi(pvs))],[std(modLow(pvs))/sqrt(length(pvs)) std(modHi(pvs))/sqrt(length(pvs))],'LineWidth',2)
xlim([0.8 2.2])
ylim([-1 1])
xlabel('3 mW vs 15 mW')
ylabel('Modulation Index')
set(gca,'TickDir','out');
title(strcat('msns Mod Index. Pval:',num2str(pvalpv)))

spikeGraphName = 'pvAndmsnsModCompare';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
