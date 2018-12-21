%This is code to plot out walls of spikes. 

cd Z:\Max\ProjectFSIAudStr\OriginalPVNpHR\PVNphRAnalysisFiles

load('180918_ML180618E_R_AudStr_3300_3mWPVHaloTuningWhiteAltLaserFullTuningAnalysis.mat', 's')

testData = s.nt10cluster3;

% hFig = figure;
spacerLength = 10;
histStore = [];

lengthTrace = length([80:120]);
for i = 1:3
    counter = 1;
    for j = 2:17
        histStore(counter:counter + lengthTrace - 1,i) = testData.FreqDBHistograms(j,i,[80:120]);
        counter = counter + lengthTrace;
        histStore(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
    
end

histStoreLaser = [];

for i = 1:3
    counter = 1;
    for j = 2:17
        histStoreLaser(counter:counter + lengthTrace - 1,i) = testData.FreqDBHistogramsLaser(j,i,[80:120]);
        counter = counter + lengthTrace;
        histStoreLaser(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
    
end

maxVal = max(max(max(histStore)),max(max(histStoreLaser)));

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
hold on
for i = 1:3
    plot(histStore(:,i)/maxVal + i-1,'k')
    plot(histStoreLaser(:,i)/maxVal+ i-1,'g')
end

spikeGraphName = 'FSIhistogramWall';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets also plot out the overall histogram
hFig = figure;
set(hFig,'Position', [10 80 1240 850])
hold on
plot(testData.HistBinVector(40:140),squeeze(mean(testData.FreqDBHistograms(2:16,3,40:140))),'b')
plot(testData.HistBinVector(40:140),squeeze(mean(testData.FreqDBHistogramsLaser(2:16,3,40:140))),'g')
xlim([-0.2 0.3])
spikeGraphName = 'FSIOverallHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



testData = s.nt11cluster1;

% hFig = figure;
spacerLength = 10;
histStore = [];

lengthTrace = length([80:120]);
for i = 1:3
    counter = 1;
    for j = 2:17
        histStore(counter:counter + lengthTrace - 1,i) = testData.FreqDBHistograms(j,i,[80:120]);
        counter = counter + lengthTrace;
        histStore(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
    
end

histStoreLaser = [];

for i = 1:3
    counter = 1;
    for j = 2:17
        histStoreLaser(counter:counter + lengthTrace - 1,i) = testData.FreqDBHistogramsLaser(j,i,[80:120]);
        counter = counter + lengthTrace;
        histStoreLaser(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
    
end

maxVal = max(max(max(histStore)),max(max(histStoreLaser)));

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
hold on
for i = 1:3
    plot(histStore(:,i)/maxVal + i-1,'k')
    plot(histStoreLaser(:,i)/maxVal+ i-1,'g')
end

spikeGraphName = 'MSNhistogramWall';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



testData = s.nt13cluster1;

% hFig = figure;
spacerLength = 10;
histStore = [];

lengthTrace = length([80:120]);
for i = 1:3
    counter = 1;
    for j = 2:17
        histStore(counter:counter + lengthTrace - 1,i) = testData.FreqDBHistograms(j,i,[80:120]);
        counter = counter + lengthTrace;
        histStore(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
    
end

histStoreLaser = [];

for i = 1:3
    counter = 1;
    for j = 2:17
        histStoreLaser(counter:counter + lengthTrace - 1,i) = testData.FreqDBHistogramsLaser(j,i,[80:120]);
        counter = counter + lengthTrace;
        histStoreLaser(counter:counter+10-1,i) = NaN;
        counter = counter + 10;
    end
    
end

maxVal = max(max(max(histStore)),max(max(histStoreLaser)));

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
hold on
for i = 1:3
    plot(histStore(:,i)/maxVal + i-1,'k')
    plot(histStoreLaser(:,i)/maxVal+ i-1,'g')
end

spikeGraphName = 'MSN2histogramWall';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets also plot out the overall histogram
hFig = figure;
set(hFig,'Position', [10 80 1240 850])
hold on
plot(testData.HistBinVector(40:140),squeeze(mean(testData.FreqDBHistograms(2:16,3,40:140))),'b')
plot(testData.HistBinVector(40:140),squeeze(mean(testData.FreqDBHistogramsLaser(2:16,3,40:140))),'g')
xlim([-0.2 0.3])
spikeGraphName = 'MSN2OverallHist';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')