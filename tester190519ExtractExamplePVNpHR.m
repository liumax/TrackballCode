%This code is meant to extract examples from 180918_ML180618E_R_AudStr_3300_3mWPVHaloTuningWhiteAltLaserFullTuningAnalysis


load('180918_ML180618E_R_AudStr_3300_3mWPVHaloTuningWhiteAltLaserFullTuningAnalysis.mat')


%FSI is 10 cluster 3. This is a partial suppression.

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt10cluster3.AllHistograms,'k','LineWidth',2)
plot(s.nt10cluster3.HistBinVector,s.nt10cluster3.AllHistogramsLaser,'g','LineWidth',2)

set(gca,'TickDir','out')
spikeGraphName = 'ExampleFSISuppression';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt10cluster3.AllHistograms,'k','LineWidth',2)
% plot(s.nt10cluster3.HistBinVector,s.nt10cluster3.AllHistogramsLaser,'g')

set(gca,'TickDir','out')
spikeGraphName = 'ExampleFSIPreSuppression';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
hold on
tarDB = 2;
maxVal = max([max(s.nt10cluster3.BinTone(2:end,tarDB)),max(s.nt10cluster3.BinToneLaser(2:end,tarDB))]);
plot(s.nt10cluster3.BinTone(2:end,tarDB),'k','LineWidth',2)
plot(s.nt10cluster3.BinToneLaser(2:end,tarDB),'g','LineWidth',2)
set(gca,'TickDir','out')
spikeGraphName = 'ExampleFSITuningCurve';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
hold on
tarDB = 2;
maxVal = max([max(s.nt10cluster3.BinTone(2:end,tarDB)),max(s.nt10cluster3.BinToneLaser(2:end,tarDB))]);
plot(s.nt10cluster3.BinTone(2:end,tarDB),s.nt10cluster3.BinToneLaser(2:end,tarDB),'k.')
[b,bintr,bintjm] = gmregress(s.nt10cluster3.BinTone(2:end,tarDB),s.nt10cluster3.BinToneLaser(2:end,tarDB));
plot(s.nt10cluster3.BinTone(2:end,tarDB),s.nt10cluster3.BinTone(2:end,tarDB)*b(2) + b(1),'r')
plot([0 maxVal],[0 maxVal],'k')
axis square
xlim([0 maxVal])
ylim([0 maxVal])
set(gca,'TickDir','out')
spikeGraphName = 'ExampleFSIScatterPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')





%plot MSNs?
hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt11cluster1.AllHistograms,'k')
plot(s.nt10cluster3.HistBinVector,s.nt11cluster1.AllHistogramsLaser,'g')

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt11cluster3.AllHistograms,'k')
plot(s.nt10cluster3.HistBinVector,s.nt11cluster3.AllHistogramsLaser,'g')

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt13cluster1.AllHistograms,'k','LineWidth',2)
plot(s.nt10cluster3.HistBinVector,s.nt13cluster1.AllHistogramsLaser,'g','LineWidth',2)
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNScalingHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt13cluster1.AllHistograms,'k','LineWidth',2)
% plot(s.nt10cluster3.HistBinVector,s.nt13cluster1.AllHistogramsLaser,'g')
ylim([0 25])
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNScalingPreHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Generate wall of histograms for selected MSN

testData = s.nt13cluster1;
histWall = [];

for i = 1:3
    counter = 1;
    for j = 2:17
        histWall(i,counter:counter + 30) = smooth(squeeze(testData.FreqDBHistograms(j,i,81:111)),5);
        counter = counter + 31;
        histWall(i,counter: counter + 9) = NaN(10,1);
        counter = counter + 10;
    end
end

histWallLaser = [];

for i = 1:3
    counter = 1;
    for j = 2:17
        histWallLaser(i,counter:counter + 30) = smooth(squeeze(testData.FreqDBHistogramsLaser(j,i,81:111)),5);
        counter = counter + 31;
        histWallLaser(i,counter: counter + 9) = NaN(10,1);
        counter = counter + 10;
    end
end

truemax = max([max(max(histWall)),max(max(histWallLaser))]);

hFig = figure;
subplot(2,1,1)
hold on
for i = 1:3
    plot(histWall(i,:)/truemax + 1*(i-1),'k')
end
set(gca,'TickDir','out')
title(num2str(truemax))
subplot(2,1,2)
hold on
for i = 1:3
    plot(histWallLaser(i,:)/truemax + 1*(i-1),'k')
end
set(gca,'TickDir','out')


spikeGraphName = 'ExampleMSNWallOfHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
plot(testData.BinDiff(2:end,2,2),'k')
hold on
plot(testData.BinDiffLaser(2:end,2,2),'g')
spikeGraphName = 'ExampleMSNTuningCurveOverlay';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')




%%




hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt13cluster2.AllHistograms,'k','LineWidth',2)
plot(s.nt10cluster3.HistBinVector,s.nt13cluster2.AllHistogramsLaser,'g','LineWidth',2)
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNAdditionHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt13cluster2.AllHistograms,'k','LineWidth',2)
% plot(s.nt10cluster3.HistBinVector,s.nt13cluster2.AllHistogramsLaser,'g')
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNAdditionPreHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0') 

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt13cluster3.AllHistograms,'k')
plot(s.nt10cluster3.HistBinVector,s.nt13cluster3.AllHistogramsLaser,'g')

hFig = figure;
hold on
plot(s.nt10cluster3.HistBinVector,s.nt14cluster1.AllHistograms,'k')
plot(s.nt10cluster3.HistBinVector,s.nt14cluster1.AllHistogramsLaser,'g')


hFig = figure;
hold on
tarDB = 2;
maxVal = max([max(s.nt13cluster1.BinTone(2:end,tarDB)),max(s.nt13cluster1.BinToneLaser(2:end,tarDB))]);
plot(s.nt13cluster1.BinTone(2:end,tarDB),'k','LineWidth',2)
plot(s.nt13cluster1.BinToneLaser(2:end,tarDB),'g','LineWidth',2)
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNScalingTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
hold on
tarDB = 2;
maxVal = max([max(s.nt13cluster1.BinTone(2:end,tarDB)),max(s.nt13cluster1.BinToneLaser(2:end,tarDB))]);
plot(s.nt13cluster1.BinTone(2:end,tarDB),s.nt13cluster1.BinToneLaser(2:end,tarDB),'k.')
[b,bintr,bintjm] = gmregress(s.nt13cluster1.BinTone(2:end,tarDB),s.nt13cluster1.BinToneLaser(2:end,tarDB));
plot(s.nt13cluster1.BinTone(2:end,tarDB),s.nt13cluster1.BinTone(2:end,tarDB)*b(2) + b(1),'r')
plot([0 maxVal],[0 maxVal],'k')
axis square
xlim([0 maxVal])
ylim([0 maxVal])
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNScalingScatter';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
hold on
tarDB = 2;
maxVal = max([max(s.nt13cluster2.BinTone(2:end,tarDB)),max(s.nt13cluster2.BinToneLaser(2:end,tarDB))]);
plot(s.nt13cluster2.BinTone(2:end,tarDB),'k','LineWidth',2)
plot(s.nt13cluster2.BinToneLaser(2:end,tarDB),'g','LineWidth',2)
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNAdditionTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



hFig = figure;
hold on
tarDB = 2;
maxVal = max([max(s.nt13cluster2.BinTone(2:end,tarDB)),max(s.nt13cluster2.BinToneLaser(2:end,tarDB))]);
plot(s.nt13cluster2.BinTone(2:end,tarDB),s.nt13cluster2.BinToneLaser(2:end,tarDB),'k.')
[b,bintr,bintjm] = gmregress(s.nt13cluster2.BinTone(2:end,tarDB),s.nt13cluster2.BinToneLaser(2:end,tarDB));
plot(s.nt13cluster2.BinTone(2:end,tarDB),s.nt13cluster2.BinTone(2:end,tarDB)*b(2) + b(1),'r')
plot([0 maxVal],[0 maxVal],'k')
axis square
xlim([0 maxVal])
ylim([0 maxVal])
set(gca,'TickDir','out')
spikeGraphName = 'ExampleMSNAdditionScatter';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

