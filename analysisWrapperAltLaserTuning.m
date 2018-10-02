% analysisWrapperAltLaserTuning




%identify files, pull names, set up for loop for extraction

% clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);

%generate empty things here/store variables

%variables
sigThresh = 5; %minimum number of significant responses at 0.05 p threshold. 
tarWin = 2; %target window of analysis. 1 = fast, 2 = tone, 3 = gen. 
%counters of cells overall
numCells = 0;
numMSNs = 0;
numPVs = 0;

%counter of number of cells I'm keeping
fullCounter = 1;
fullNames = []; %for storage of file and unit names
fullMaster = [];
fullWidth = []; %for storage of width values
fullBin = [];
fullSig = [];
fullLat = [];



for i = 1:numFiles
    %load the target file
    load(targetFiles{i})
    numUnits = size(masterData,1);
    numDBs = s.SoundData.NumDBs;
    %update total number of cells
    numCells = numCells + numUnits;
    %pull number of MSN and PVs
    numMSNs = numMSNs + length(find(masterData(:,7) == 0));
    numPVs = numPVs + length(find(masterData(:,7) == 1));
    %now I want to pick out only significantly responsive cells. 
    %here, I pull out the data from binSigVals, since the width data is for
    %contiguous segments, and therefore ignores many units. 
    keepList = zeros(numUnits,1);
    for j = 1:numUnits
        for k = 1:numDBs
            sigVals(k) = length(find(s.(s.DesignationName{j}).BinSigVals(:,:,k) < 0.05));
            sigValsLaser(k) = length(find(s.(s.DesignationName{j}).BinSigValsLaser(:,:,k) < 0.05));
        end
        maxVal = max([max(sigVals),max(sigValsLaser)]);
        if maxVal > sigThresh
            keepList(j) = 1;
        end
    end
    %now that I've found significantly responding units, lets keep just
    %these. I'll need to extract relevant information about these units,
    %including their identity (PV vs MSN vs UNK), tuning curve width (with
    %and without laser) at all windows, binned spikes (with and without
    %laser) at all windows, using mean subtracted, latency times (based on
    %original latency code). Also keep mean firing rate, peak trough time. 
    keepList = find(keepList == 1);
    for j = 1:length(keepList)
        %store name
        trueName = strcat(targetFiles{i}(1:end-4),s.DesignationName{keepList(j)});
        fullNames{fullCounter} = trueName;
        %first store designation and basic info
        fullMaster(fullCounter,1) = masterData(keepList(j),7); %store PV/MSN desig
        fullMaster(fullCounter,2) = masterData(keepList(j),5); %store peak trough time
        fullMaster(fullCounter,3) = masterData(keepList(j),4); %store average rate
        fullMaster(fullCounter,4) = mean(diff(s.SoundData.UniqueDBs));%store db step size
        %now pull width values for positive responses
        fullWidth(fullCounter,:,:) = [squeeze(s.NonLaserOverall.PosWidths(:,keepList(j),:)),squeeze(s.LaserOverall.PosWidths(:,keepList(j),:))]; 
        %stores a 3 x 6 matrix. columns are within a window (going from fast, tone, gen) with first three being no-laser, last three being laser. rows are DB increases
        
        %store which values in target window have significant response. 
        findSigTone = find(s.(s.DesignationName{keepList(j)}).BinSigVals(:,:,tarWin)<0.05);
        testMask = zeros(size(s.(s.DesignationName{keepList(j)}).BinSigVals(:,:,tarWin)));
        testMask(findSigTone) = 1;
        fullSig{fullCounter} = testMask;
        
        %now pull binned spike values from both laser and non-laser
        fullBin{fullCounter} = [s.(s.DesignationName{keepList(j)}).BinDiff(:,:,tarWin),s.(s.DesignationName{keepList(j)}).BinDiffLaser(:,:,tarWin)];
        %now do some calculations
        binDiffs = (s.(s.DesignationName{keepList(j)}).BinDiffLaser(:,:,tarWin) - s.(s.DesignationName{keepList(j)}).BinDiff(:,:,tarWin)).*testMask;
        %also try and pull a modulation index of response
        toneMod = (nanmean(s.(s.DesignationName{keepList(j)}).BinDiffLaser(:,:,tarWin).*testMask) - nanmean(s.(s.DesignationName{keepList(j)}).BinDiff(:,:,tarWin).*testMask))/(nanmean(s.(s.DesignationName{keepList(j)}).BinDiffLaser(:,:,tarWin).*testMask) + nanmean(s.(s.DesignationName{keepList(j)}).BinDiff(:,:,tarWin).*testMask));
        toneModStore(fullCounter) = toneMod;
        realDiffs = binDiffs(binDiffs ~= 0);
        fullMaster(fullCounter,5) = nanmean(realDiffs);
        fullMaster(fullCounter,6) = length(realDiffs);
        curveStore{fullCounter} = s.(s.DesignationName{keepList(j)}).BinDiff(:,:,tarWin);
        curveStoreLaser{fullCounter} = s.(s.DesignationName{keepList(j)}).BinDiffLaser(:,:,tarWin);
        %now store latency data. 
        fullLat{fullCounter} = [s.(s.DesignationName{keepList(j)}).LatencyMap,s.(s.DesignationName{keepList(j)}).LatencyMapLaser];
        latDiffs = (s.(s.DesignationName{keepList(j)}).LatencyMapLaser - s.(s.DesignationName{keepList(j)}).LatencyMap).*testMask;
        realDiffs = latDiffs(latDiffs ~= 0);
        fullMaster(fullCounter,7) = nanmean(realDiffs);
        fullMaster(fullCounter,8) = length(realDiffs);
        fullMaster(fullCounter,9) = length(findSigTone);
        
        %now store data from positive widths
        fullMaster(fullCounter,10) = s.NonLaserOverall.PosWidths(3,keepList(j),tarWin);
        fullMaster(fullCounter,11) = s.LaserOverall.PosWidths(3,keepList(j),tarWin);
        
        %store modulation indices
        fullMaster(fullCounter,12) = masterData(keepList(j),22);
        
        %now store positive width data from lower amplitudes
        fullMaster(fullCounter,13) = s.NonLaserOverall.PosWidths(1,keepList(j),tarWin);
        fullMaster(fullCounter,14) = s.LaserOverall.PosWidths(1,keepList(j),tarWin);
        fullMaster(fullCounter,15) = s.NonLaserOverall.PosWidths(2,keepList(j),tarWin);
        fullMaster(fullCounter,16) = s.LaserOverall.PosWidths(2,keepList(j),tarWin);
        
        fullCounter = fullCounter + 1;
    end
end




msns = find(fullMaster(:,1) == 0);
pvs = find(fullMaster(:,1) == 1);



nanmean(fullMaster(msns,5))
nanmean(fullMaster(msns,7))

nanmean(fullMaster(pvs,5))
nanmean(fullMaster(pvs,7))


hFig = figure;
set(hFig, 'Position', [10 80 1000 1000])
subplot(2,3,1)
hist(fullMaster(msns,12),[-1:0.1:1])
xlim([-1 1])
xlabel('Modulation Index')
ylabel('Number MSNs')
title('Modulation Index MSNs')

subplot(2,3,4)
hist(fullMaster(pvs,12),[-1:0.1:1])
xlim([-1 1])
xlabel('Modulation Index')
ylabel('Number PVs')
title('Modulation Index PVs')

subplot(2,3,2)
hist(toneModStore(msns),[-1:0.1:1])
xlim([-1 1])
xlabel('Modulation Index')
ylabel('Number MSNs')
title('Tone Modulation Index MSNs')

subplot(2,3,5)
hist(toneModStore(pvs),[-1:0.1:1])
xlim([-1 1])
xlabel('Modulation Index')
ylabel('Number PVs')
title('Tone Modulation Index PVs')

subplot(2,3,3)
plot(fullMaster(msns,12),toneModStore(msns),'k.')
hold on
plot([-1 1],[-1 1],'r')
xlabel('ModulationIndexBase')
ylabel('ModulationIndexTone')
title('MSN ModIndex Comparison')

subplot(2,3,6)
plot(fullMaster(pvs,12),toneModStore(pvs),'k.')
hold on
plot([-1 1],[-1 1],'r')
xlabel('ModulationIndexBase')
ylabel('ModulationIndexTone')
title('PV ModIndex Comparison')


spikeGraphName = 'Modulation Indices';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
set(hFig, 'Position', [10 80 1000 1000])
subplot(2,1,1)
hold on
plot(fullMaster(msns,6),fullMaster(msns,5),'k.')
plot([0 max(fullMaster(msns,6))],[0 0],'r')
xlabel('Number of Non-Zero Responses')
ylabel('Mean Change in Spikes/Tone')
title('Effect of Laser on Size of Responses')

subplot(2,1,2)
hold on
plot(fullMaster(msns,13),fullMaster(msns,14),'k.')
plot(fullMaster(msns,15),fullMaster(msns,16),'bo')
plot(fullMaster(msns,10),fullMaster(msns,11),'m*')
plot([0 max(max(fullMaster(msns,10:11)))],[0 max(max(fullMaster(msns,10:11)))],'b')
xlabel('Tuning Width No Laser')
ylabel('Tuning Width With Laser')
title('Effect of Laser on Width of Tuning Curve at 70 dB SPL')

spikeGraphName = 'PV NpHR Effects';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at the average change in width at different amplitude levels. I
%want to look at this both as a ratio and as an absolute difference
ratioLow = fullMaster(msns,14)./fullMaster(msns,13);
ratioMed = fullMaster(msns,16)./fullMaster(msns,15);
ratioHi = fullMaster(msns,11)./fullMaster(msns,10);

ratioLow(ratioLow == Inf) = NaN;
ratioMed(ratioMed == Inf) = NaN;
ratioHi(ratioHi == Inf) = NaN;

widthDiff(:,1)= fullMaster(msns,14)-fullMaster(msns,13);
widthDiff(:,2)= fullMaster(msns,16)-fullMaster(msns,15);
widthDiff(:,3)= fullMaster(msns,11)-fullMaster(msns,10);

%plot out subtraction values
hFig = figure;
subplot(3,2,1)
hist(widthDiff(:,1),[-5:1:10])
xlim([-5 10])
subplot(3,2,3)
hist(widthDiff(:,2),[-5:1:10])
xlim([-5 10])
subplot(3,2,5)
hist(widthDiff(:,3),[-5:1:10])
xlim([-5 10])
subplot(1,2,2)
plot(widthDiff')


%plot out tuning curves for MSNs
hFig = figure;
ampLevel = 1;
set(hFig, 'Position', [10 80 1000 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(msns)
    subplot(10,10,i)
    hold on
    plot([1 length(curveStore{msns(i)})-1],[0 0],'k')
    plot(curveStore{msns(i)}(2:end,ampLevel),'r')
    plot(curveStoreLaser{msns(i)}(2:end,ampLevel),'b')
    xlim([1,length(curveStore{msns(i)})-1])
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
end

hFig = figure;
% ampLevel = 3;
set(hFig, 'Position', [10 80 1000 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(pvs)
    subplot(6,6,i)
    hold on
    plot([1 length(curveStore{pvs(i)})-1],[0 0],'k')
    plot(curveStore{pvs(i)}(2:end,ampLevel),'r')
    plot(curveStoreLaser{pvs(i)}(2:end,ampLevel),'b')
    xlim([1,length(curveStore{pvs(i)})-1])
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
end

%plot heatmap tuning curves
hFig = figure;
% ampLevel = 2;
set(hFig, 'Position', [10 80 1000 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(msns)
    subplot(10,10,i)
    imagesc(curveStore{msns(i)}(2:end,:)')
    colormap('parula')
%     hold on
%     plot([1 length(curveStore{msns(i)})-1],[0 0],'k')
%     plot(curveStore{msns(i)}(2:end,ampLevel),'r')
%     plot(curveStoreLaser{msns(i)}(2:end,ampLevel),'b')
%     xlim([1,length(curveStore{msns(i)})-1])
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
end

hFig = figure;
% ampLevel = 2;
set(hFig, 'Position', [10 80 1000 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(pvs)
    subplot(6,6,i)
    imagesc(curveStore{pvs(i)}(2:end,:)')
    colormap('parula')
%     hold on
%     plot([1 length(curveStore{msns(i)})-1],[0 0],'k')
%     plot(curveStore{msns(i)}(2:end,ampLevel),'r')
%     plot(curveStoreLaser{msns(i)}(2:end,ampLevel),'b')
%     xlim([1,length(curveStore{msns(i)})-1])
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
end







