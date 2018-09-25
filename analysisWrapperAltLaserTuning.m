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
sigThresh = 3; %minimum number of significant responses at 0.05 p threshold. 
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
        realDiffs = binDiffs(binDiffs ~= 0);
        fullMaster(fullCounter,5) = nanmean(realDiffs);
        fullMaster(fullCounter,6) = length(realDiffs);
        %now store latency data. 
        fullLat{fullCounter} = [s.(s.DesignationName{keepList(j)}).LatencyMap,s.(s.DesignationName{keepList(j)}).LatencyMapLaser];
        latDiffs = (s.(s.DesignationName{keepList(j)}).LatencyMapLaser - s.(s.DesignationName{keepList(j)}).LatencyMap).*testMask;
        realDiffs = latDiffs(latDiffs ~= 0);
        fullMaster(fullCounter,7) = nanmean(realDiffs);
        fullMaster(fullCounter,8) = length(realDiffs);
        fullMaster(fullCounter,9) = length(findSigTone);
        
        %now store data from positive widths
        fullMaster(fullCounter,10) = s.NonLaserOverall.PosWidths(3,keepList(j),3);
        fullMaster(fullCounter,11) = s.LaserOverall.PosWidths(3,keepList(j),3);
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
subplot(2,1,1)
hold on
plot(fullMaster(msns,6),fullMaster(msns,5),'k.')
plot([0 max(fullMaster(msns,6))],[0 0],'r')
xlabel('Number of Significant Responses')
ylabel('Mean Change in Spikes/Tone')
title('Effect of Laser on Size of Responses')

subplot(2,1,2)
hold on
plot(fullMaster(msns,10),fullMaster(msns,11),'k.')
plot([0 max(max(fullMaster(msns,10:11)))],[0 max(max(fullMaster(msns,10:11)))],'r')
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














