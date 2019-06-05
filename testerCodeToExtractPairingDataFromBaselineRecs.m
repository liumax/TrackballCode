
clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'Analysis');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);


%parameters
rpvVect = [0:0.0005:0.025];
rpvLim = 0.5; %in %. 
rpvCut = 0.0015; %cutoff for considering as RPV. 
prctileBounds = [0.5 99.5];
posLims = [15:25]; %target zone for detecting positives
negLimsL = [5:18]; %target zone for detecting negatives
negLimsR = [22:35]; %target zone for detecting negatives
diffLim = 2; %limit for consecutive nature of negative deviations. 

bigMaster = [];
bigMasterInd = 1;
respVect = [];
pkTroughRatioStore = [];
binValBigStore = [];
sigValBigStore = [];
latMapBigStore = [];
widthLatStore = [];
nameStore = [];
unitStore = [];
halfWidthTimeStore = [];
pkTroughTimeStore = [];

%create store for types of interactions in cross correlogram. 
%ones I can think of are: Positive (+/-), and Negative (+ pre, + post,+ both,  -).
%This makes 4 categories. 

intTypeStore = zeros(4,2,3); %dim1 = negative (negative, + pre, + post, + both), dim2 = positive (neg, +), and dim 3 = type (FSI -> MSN, MSN -> MSN FSI -> FSI)

%parameters
masterHeaderSize = 11; %only want the first 12 values of masterData. 


bigInfo = zeros(1,10);
bigxCorr = [];
bigbaseCorr = [];
bigCounter = 1;
bigSigNeg = [];
bigSigPos = [];

%% Extract correlation information from files. 
for i = 1:numFiles
    %load the target file
    disp(strcat('Analyzing File:',num2str(i)))
    load(targetFiles{i})
    
    %determine number of units
    numUnits = size(masterData,1);
    %determine unit names
    cluNames = s.DesignationName;
    
    %now march through each shank 
    shank1Units = find(s.PeakChanVals <=32);
    subUnits = length(shank1Units);
    counter = 1;
    infoStoreShank1 = zeros(1,8);
    for j = 1:subUnits
        %venture through the units.Go up until the unit before the same unit.  
        for k = 1:j-1
            %store information about aligning unit
             infoStoreShank1(counter,1) = shank1Units(j);
             infoStoreShank1(counter,2) = masterData(shank1Units(j),6);
             %store information about aligned unit
             infoStoreShank1(counter,3) = shank1Units(k);
             infoStoreShank1(counter,4) = masterData(shank1Units(k),6);
             %store distance
             infoStoreShank1(counter,5) = abs(s.PeakChanVals(shank1Units(j))-s.PeakChanVals(shank1Units(k)))*25;
             %now we want to look for significance. First generate the
             %percentile bounds
             tmpPercentiles = prctile(corrData.ShuffStore1{j,k},prctileBounds,2);
             tmpTrueVal = squeeze(corrData.trueStoreShank1(j,k,:));
             posVals = tmpTrueVal - tmpPercentiles(:,2);
             negVals = tmpTrueVal - tmpPercentiles(:,1);
             prctileStore(counter,:,:) = tmpTrueVal - tmpPercentiles;
             infoStoreShank1(counter,6) = length(find(posVals(posLims) > 0));
             infoStoreShank1(counter,7) = length(find(negVals(negLimsL) < 0));
             infoStoreShank1(counter,8) = length(find(negVals(negLimsR) < 0));
             xcorrStore(counter,:) = tmpTrueVal;
             basecorrStore(counter,:) = mean(corrData.ShuffStore1{j,k}');
             counter = counter + 1;
        end
    end
    
    %put into big store!
    bigInfo(bigCounter:bigCounter + size(infoStoreShank1,1) - 1,1) = i;
    bigInfo(bigCounter:bigCounter + size(infoStoreShank1,1) - 1,2) = 1;
    bigInfo(bigCounter:bigCounter + size(infoStoreShank1,1) - 1,3:10) = infoStoreShank1;
    bigxCorr(bigCounter:bigCounter + size(infoStoreShank1,1) - 1,:) = xcorrStore;
    bigbaseCorr(bigCounter:bigCounter + size(infoStoreShank1,1) - 1,:) = basecorrStore;
    bigSigNeg(bigCounter:bigCounter + size(infoStoreShank1,1) - 1,:) = prctileStore(:,:,1);
    bigSigPos(bigCounter:bigCounter + size(infoStoreShank1,1) - 1,:) = prctileStore(:,:,2);
    bigCounter = bigCounter + size(infoStoreShank1,1);
    xcorrStore = [];
    basecorrStore = [];
    prctileStore = [];
    
    %shank 2
    shank2Units = find(s.PeakChanVals > 32);
    subUnits = length(shank2Units);
    counter = 1;
    infoStoreShank2 = zeros(1,8);
    for j = 1:subUnits
        %venture through the units.Go up until the unit before the same unit.  
        for k = 1:j-1
            %store information about aligning unit
             infoStoreShank2(counter,1) = shank2Units(j);
             infoStoreShank2(counter,2) = masterData(shank2Units(j),6);
             %store information about aligned unit
             infoStoreShank2(counter,3) = shank2Units(k);
             infoStoreShank2(counter,4) = masterData(shank2Units(k),6);
             %store distance
             infoStoreShank2(counter,5) = abs(s.PeakChanVals(shank2Units(j))-s.PeakChanVals(shank2Units(k)))*25;
             %now we want to look for significance. First generate the
             %percentile bounds
             tmpPercentiles = prctile(corrData.ShuffStore2{j,k},prctileBounds,2);
             tmpTrueVal = squeeze(corrData.trueStoreShank2(j,k,:));
             posVals = tmpTrueVal - tmpPercentiles(:,2);
             negVals = tmpTrueVal - tmpPercentiles(:,1);
             prctileStore(counter,:,:) = tmpTrueVal - tmpPercentiles;
             infoStoreShank2(counter,6) = length(find(posVals(posLims) > 0));
             infoStoreShank2(counter,7) = length(find(negVals(negLimsL) < 0));
             infoStoreShank2(counter,8) = length(find(negVals(negLimsR) < 0));
             xcorrStore(counter,:) = tmpTrueVal;
             basecorrStore(counter,:) = mean(corrData.ShuffStore2{j,k}');
             counter = counter + 1;
        end
    end
    
    %put into big store!
    bigInfo(bigCounter:bigCounter + size(infoStoreShank2,1) - 1,1) = i;
    bigInfo(bigCounter:bigCounter + size(infoStoreShank2,1) - 1,2) = 2;
    bigInfo(bigCounter:bigCounter + size(infoStoreShank2,1) - 1,3:10) = infoStoreShank2;
    bigxCorr(bigCounter:bigCounter + size(infoStoreShank2,1) - 1,:) = xcorrStore;
    bigbaseCorr(bigCounter:bigCounter + size(infoStoreShank2,1) - 1,:) = basecorrStore;
    bigSigNeg(bigCounter:bigCounter + size(infoStoreShank2,1) - 1,:) = prctileStore(:,:,1);
    bigSigPos(bigCounter:bigCounter + size(infoStoreShank2,1) - 1,:) = prctileStore(:,:,2);
    bigCounter = bigCounter + size(infoStoreShank2,1);
    xcorrStore = [];
    basecorrStore = [];
    prctileStore = [];
    
end

%% Extract waveform information from files
%we also want to get parameters of waveforms
%actually extract files.
waveInfoStore = zeros(1,2);
waveFormStore = [];
waveCount = 1;
for i = 1:numFiles
    %load the target file
    disp(strcat('Analyzing File:',num2str(i)))
    load(targetFiles{i})
    %determine number of units
    numUnits = size(masterData,1);
    %determine unit names
    cluNames = s.DesignationName;
    
    waveInfoStore(waveCount:waveCount + numUnits - 1,:) = masterData(:,1:2);
    %go through and save waveforms. 
    for j = 1:numUnits
        waveFormStore(:,:,waveCount + j - 1) = s.(cluNames{j}).medianWave;
    end
    
    waveCount = waveCount + numUnits;
    
end


%% Process correlation information: remove NaN cell types
%now we want to eliminate unclear cell types. 
findNan1 = find(isnan(bigInfo(:,4)));
findNan2 = find(isnan(bigInfo(:,6)));
uniqueNan = unique([findNan1;findNan2]);
bigInfoBackup = bigInfo; %stash a backup. 

%delete NaN from bigInfo
bigInfo(uniqueNan,:) = NaN;
%now go through and pull only the units within 250
findFar = find(bigInfo(:,7) > 250);
bigInfo(findFar,:) = NaN;

%find where its 1 and then 0 for cell types
cellTypeFind = bigInfo(:,4) - bigInfo(:,6);

%% Process correlation information: Look at PV -> MSN inhibition interactions
fsiTar = 1; %tells where to store in array!
%find where aligned to PV. 
alignedPV = find(cellTypeFind == 1);
%pull the pos, neg L neg R
alignedPVSig = bigInfo(alignedPV,8:10);

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedPVSig(:,2) == 0);
findPost = find(alignedPVSig(:,3) > 1);
combFind = intersect(findPre,findPost);

%now i need to determine if negatives are consecutive or separated. 
qcCheck = zeros(length(combFind),1);
for i = 1:length(combFind)
    findNeg = find(bigSigNeg(alignedPV(combFind(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFind = combFind(qcCheck==1);

%determine number positive
findPos = find(alignedPVSig(combFind,1) >1);
%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPV(combFind(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

%store number
intTypeStore(3,1,fsiTar) = intTypeStore(3,1,fsiTar) + (length(combFind) - length(findPos));
intTypeStore(3,2,fsiTar) = intTypeStore(3,2,fsiTar) + length(findPos);

%find where reverse aligned
alignedPVRev = find(cellTypeFind == -1);
alignedPVSigRev = bigInfo(alignedPVRev,[8,10,9]); %fix order here. 

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedPVSigRev(:,2) == 0);
findPost = find(alignedPVSigRev(:,3) > 1);
combFindRev = intersect(findPre,findPost);

%now i need to determine if negatives are consecutive or separated. 
qcCheck = zeros(length(combFindRev),1);
for i = 1:length(combFindRev)
    findNeg = find(bigSigNeg(alignedPVRev(combFindRev(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFindRev = combFindRev(qcCheck == 1);
%determine number positive
findPos = find(alignedPVSigRev(combFindRev,1) >1);
%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPVRev(combFindRev(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);
%store number
intTypeStore(3,1,fsiTar) = intTypeStore(3,1,fsiTar) + (length(combFindRev) - length(findPos));
intTypeStore(3,2,fsiTar) = intTypeStore(3,2,fsiTar) + length(findPos);

% %plot out cross correlograms of putative significant ones
% 
% figure
% subPlotLim = ceil(sqrt(length(combFind)));
% for i = 1:length(combFind)
%     subplot(subPlotLim,subPlotLim,i)
%     plot(bigxCorr(alignedPV(combFind(i)),:),'k')
%     hold on
%     findNeg = find(bigSigNeg(alignedPV(combFind(i)),:) < 0);
%     plot(findNeg,bigxCorr(alignedPV(combFind(i)),findNeg),'g*')
%     findPos = find(bigSigPos(alignedPV(combFind(i)),:) > 0);
%     plot(findPos,bigxCorr(alignedPV(combFind(i)),findPos),'m*')
% end
% 
% %plot with shuffle subtracted
% figure
% subPlotLim = ceil(sqrt(length(combFind)));
% for i = 1:length(combFind)
%     subplot(subPlotLim,subPlotLim,i)
%     plot(bigxCorr(alignedPV(combFind(i)),:) - bigbaseCorr(alignedPV(combFind(i)),:),'k')
%     hold on
%     findNeg = find(bigSigNeg(alignedPV(combFind(i)),:) < 0);
%     plot(findNeg,bigxCorr(alignedPV(combFind(i)),findNeg)- bigbaseCorr(alignedPV(combFind(i)),findNeg),'g*')
%     findPos = find(bigSigPos(alignedPV(combFind(i)),:) > 0);
%     plot(findPos,bigxCorr(alignedPV(combFind(i)),findPos)- bigbaseCorr(alignedPV(combFind(i)),findPos),'m*')
% end
% 
% %plot out reverse units: cross correlograms of putative significant ones
% 
% figure
% subPlotLim = ceil(sqrt(length(combFindRev)));
% for i = 1:length(combFindRev)
%     subplot(subPlotLim,subPlotLim,i)
%     plot(bigxCorr(alignedPVRev(combFindRev(i)),:),'k')
%     hold on
%     findNeg = find(bigSigNeg(alignedPVRev(combFindRev(i)),:) < 0);
%     plot(findNeg,bigxCorr(alignedPVRev(combFindRev(i)),findNeg),'g*')
%     findPos = find(bigSigPos(alignedPVRev(combFindRev(i)),:) > 0);
%     plot(findPos,bigxCorr(alignedPVRev(combFindRev(i)),findPos),'m*')
% end
% 
% %plot with shuffle subtracted
% figure
% subPlotLim = ceil(sqrt(length(combFindRev)));
% for i = 1:length(combFindRev)
%     subplot(subPlotLim,subPlotLim,i)
%     plot(bigxCorr(alignedPVRev(combFindRev(i)),:) - bigbaseCorr(alignedPVRev(combFindRev(i)),:),'k')
%     hold on
%     findNeg = find(bigSigNeg(alignedPVRev(combFindRev(i)),:) < 0);
%     plot(findNeg,bigxCorr(alignedPVRev(combFindRev(i)),findNeg)- bigbaseCorr(alignedPVRev(combFindRev(i)),findNeg),'g*')
%     findPos = find(bigSigPos(alignedPVRev(combFindRev(i)),:) > 0);
%     plot(findPos,bigxCorr(alignedPVRev(combFindRev(i)),findPos)- bigbaseCorr(alignedPVRev(combFindRev(i)),findPos),'m*')
% end




    
%% Process correlation information: Look at PV -> MSN inhibition before zero

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedPVSig(:,3) == 0);
findPost = find(alignedPVSig(:,2) > 1);

combFindPre = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combFindPre),1);
for i = 1:length(combFindPre)
    findNeg = find(bigSigNeg(alignedPV(combFindPre(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFindPre = combFindPre(qcCheck == 1);
%determine number positive
findPos = find(alignedPVSig(combFindPre,1) >1);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPV(combFindPre(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(2,1,fsiTar) = intTypeStore(2,1,fsiTar) + (length(combFindPre) - length(findPos));
intTypeStore(2,2,fsiTar) = intTypeStore(2,2,fsiTar) + length(findPos);


%do for reversed units
findPre = find(alignedPVSigRev(:,3) == 0);
findPost = find(alignedPVSigRev(:,2) > 1);

combFindPreRev = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combFindPreRev),1);
for i = 1:length(combFindPreRev)
    findNeg = find(bigSigNeg(alignedPVRev(combFindPreRev(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFindPreRev = combFindPreRev(qcCheck == 1);
%determine number positive
findPos = find(alignedPVSigRev(combFindPreRev,1) >1);
%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPVRev(combFindPreRev(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(2,1,fsiTar) = intTypeStore(2,1,fsiTar) + (length(combFindPreRev) - length(findPos));
intTypeStore(2,2,fsiTar) = intTypeStore(2,2,fsiTar) + length(findPos);

% figure
% subPlotLim = ceil(sqrt(length(combFindPre)));
% for i = 1:length(combFindPre)
%     subplot(subPlotLim,subPlotLim,i)
%     plot(bigxCorr(alignedPV(combFindPre(i)),:),'k')
%     hold on
%     findNeg = find(bigSigNeg(alignedPV(combFindPre(i)),:) < 0);
%     plot(findNeg,bigxCorr(alignedPV(combFindPre(i)),findNeg),'g*')
%     findPos = find(bigSigPos(alignedPV(combFindPre(i)),:) > 0);
%     plot(findPos,bigxCorr(alignedPV(combFindPre(i)),findPos),'m*')
% end

%% Process correlation information: Look at PV -> MSN inhibition before and after

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedPVSig(:,3) >= 1);
findPost = find(alignedPVSig(:,2) >= 1);

combFindBoth = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combFindBoth),1);
for i = 1:length(combFindBoth)
    findNeg = find(bigSigNeg(alignedPV(combFindBoth(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFindBoth = combFindBoth(qcCheck == 1);
%determine number positive
findPos = find(alignedPVSig(combFindBoth,1) >1);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPV(combFindBoth(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(4,1,fsiTar) = intTypeStore(4,1,fsiTar) + (length(combFindBoth) - length(findPos));
intTypeStore(4,2,fsiTar) = intTypeStore(4,2,fsiTar) + length(findPos);


%do for reversed units
findPre = find(alignedPVSigRev(:,3) >= 1);
findPost = find(alignedPVSigRev(:,2) >= 1);

combFindBothRev = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combFindBothRev),1);
for i = 1:length(combFindBothRev)
    findNeg = find(bigSigNeg(alignedPVRev(combFindBothRev(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFindBothRev = combFindBothRev(qcCheck == 1);
%determine number positive
findPos = find(alignedPVSigRev(combFindBothRev,1) >1);
%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPVRev(combFindBothRev(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(4,1,fsiTar) = intTypeStore(4,1,fsiTar) + (length(combFindBothRev) - length(findPos));
intTypeStore(4,2,fsiTar) = intTypeStore(4,2,fsiTar) + length(findPos);


%% Process correlation information: Look at PV -> MSN no inhibition

findPre = find(alignedPVSig(:,3) <= 1);
findPost = find(alignedPVSig(:,2) <= 1);

combFindNone = intersect(findPre,findPost);

%determine number positive
findPos = find(alignedPVSig(combFindNone,1) >0);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPV(combFindNone(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(1,1,fsiTar) = intTypeStore(1,1,fsiTar) + (length(combFindNone) - length(findPos));
intTypeStore(1,2,fsiTar) = intTypeStore(1,2,fsiTar) + length(findPos);


%do for reversed units
findPre = find(alignedPVSigRev(:,3) <= 1);
findPost = find(alignedPVSigRev(:,2) <= 1);

combFindNoneRev = intersect(findPre,findPost);

%determine number positive
findPos = find(alignedPVSigRev(combFindNoneRev,1) >0);
%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedPVRev(combFindNoneRev(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(1,1,fsiTar) = intTypeStore(1,1,fsiTar) + (length(combFindNoneRev) - length(findPos));
intTypeStore(1,2,fsiTar) = intTypeStore(1,2,fsiTar) + length(findPos);

%% Do some plotting
figure
axisLim = ceil(sqrt(length(combFind)));
for i = 1:length(combFind)
subplot(axisLim,axisLim,i)
plot(bigxCorr(alignedPV(combFind(i)),:))
end

figure
axisLim = ceil(sqrt(length(combFind)));
for i = 1:length(combFind)
subplot(axisLim,axisLim,i)
plot(bigxCorr(alignedPV(combFind(i)),:)-bigbaseCorr(alignedPV(combFind(i)),:))
end

figure
axisLim = ceil(sqrt(length(combFindRev)));
for i = 1:length(combFindRev)
subplot(axisLim,axisLim,i)
plot(bigxCorr(alignedPV(combFindRev(i)),[end:-1:1])-bigbaseCorr(alignedPV(combFindRev(i)),[end:-1:1]))
end



figure
axisLim = ceil(sqrt(length(combFindPre)));
for i = 1:length(combFindPre)
subplot(axisLim,axisLim,i)
plot(bigxCorr(alignedPV(combFindPre(i)),:))
end

%% Now lets do the same for MSNs!
%% Process correlation information: Look at MSN -> MSN inhibition interactions
msnTar = 2; %tells where to store in array!
%find where aligned to MSN. 
alignedMSN = find(bigInfo(:,4)==0 & bigInfo(:,6) == 0);
%pull the pos, neg L neg R
alignedMSNSig = bigInfo(alignedMSN,8:10);

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedMSNSig(:,2) == 0);
findPost = find(alignedMSNSig(:,3) > 1);
combMSNFind = intersect(findPre,findPost);

%now i need to determine if negatives are consecutive or separated. 
qcCheck = zeros(length(combMSNFind),1);
for i = 1:length(combMSNFind)
    findNeg = find(bigSigNeg(alignedMSN(combMSNFind(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combMSNFind = combMSNFind(qcCheck==1);

%determine number positive
findPos = find(alignedMSNSig(combMSNFind,1) >1);
%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedMSN(combMSNFind(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

%store number
intTypeStore(3,1,msnTar) = intTypeStore(3,1,msnTar) + (length(combMSNFind) - length(findPos));
intTypeStore(3,2,msnTar) = intTypeStore(3,2,msnTar) + length(findPos);

%% Process correlation information: Look at MSN -> MSN inhibition before zero

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedMSNSig(:,3) == 0);
findPost = find(alignedMSNSig(:,2) > 1);

combMSNFindPre = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combMSNFindPre),1);
for i = 1:length(combMSNFindPre)
    findNeg = find(bigSigNeg(alignedMSN(combMSNFindPre(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combMSNFindPre = combMSNFindPre(qcCheck == 1);
%determine number positive
findPos = find(alignedMSNSig(combMSNFindPre,1) >1);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedMSN(combMSNFindPre(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(2,1,msnTar) = intTypeStore(2,1,msnTar) + (length(combMSNFindPre) - length(findPos));
intTypeStore(2,2,msnTar) = intTypeStore(2,2,msnTar) + length(findPos);


%% Process correlation information: Look at MSN -> MSN inhibition before and after

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedMSNSig(:,3) >= 1);
findPost = find(alignedMSNSig(:,2) >= 1);

combMSNFindBoth = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combMSNFindBoth),1);
for i = 1:length(combMSNFindBoth)
    findNeg = find(bigSigNeg(alignedMSN(combMSNFindBoth(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combMSNFindBoth = combMSNFindBoth(qcCheck == 1);
%determine number positive
findPos = find(alignedMSNSig(combMSNFindBoth,1) >1);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedMSN(combMSNFindBoth(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(4,1,msnTar) = intTypeStore(4,1,msnTar) + (length(combMSNFindBoth) - length(findPos));
intTypeStore(4,2,msnTar) = intTypeStore(4,2,msnTar) + length(findPos);


%% Process correlation information: Look at MSN -> MSN no inhibition

findPre = find(alignedMSNSig(:,3) <= 1);
findPost = find(alignedMSNSig(:,2) <= 1);

combMSNFindNone = intersect(findPre,findPost);

%determine number positive
findPos = find(alignedMSNSig(combMSNFindNone,1) >0);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedMSN(combMSNFindNone(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(1,1,msnTar) = intTypeStore(1,1,msnTar) + (length(combMSNFindNone) - length(findPos));
intTypeStore(1,2,msnTar) = intTypeStore(1,2,msnTar) + length(findPos);


%% Now lets see what we can pull from FSI FSI
%% Process correlation information: Look at FSFS -> FSFS inhibition interactions
fsfsTar = 3; %tells where to store in array!
%find where aligned to FSFS. 
alignedFSFS = find(bigInfo(:,4)==1 & bigInfo(:,6) == 1);
%pull the pos, neg L neg R
alignedFSFSSig = bigInfo(alignedFSFS,8:10);

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedFSFSSig(:,2) == 0);
findPost = find(alignedFSFSSig(:,3) > 1);
combFSFSFind = intersect(findPre,findPost);

%now i need to determine if negatives are consecutive or separated. 
qcCheck = zeros(length(combFSFSFind),1);
for i = 1:length(combFSFSFind)
    findNeg = find(bigSigNeg(alignedFSFS(combFSFSFind(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFSFSFind = combFSFSFind(qcCheck==1);

%determine number positive
findPos = find(alignedFSFSSig(combFSFSFind,1) >1);
%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedFSFS(combFSFSFind(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

%store number
intTypeStore(3,1,fsfsTar) = intTypeStore(3,1,fsfsTar) + (length(combFSFSFind) - length(findPos));
intTypeStore(3,2,fsfsTar) = intTypeStore(3,2,fsfsTar) + length(findPos);

%% Process correlation information: Look at FSFS -> FSFS inhibition before zero

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedFSFSSig(:,3) == 0);
findPost = find(alignedFSFSSig(:,2) > 1);

combFSFSFindPre = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combFSFSFindPre),1);
for i = 1:length(combFSFSFindPre)
    findNeg = find(bigSigNeg(alignedFSFS(combFSFSFindPre(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFSFSFindPre = combFSFSFindPre(qcCheck == 1);
%determine number positive
findPos = find(alignedFSFSSig(combFSFSFindPre,1) >1);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedFSFS(combFSFSFindPre(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(2,1,fsfsTar) = intTypeStore(2,1,fsfsTar) + (length(combFSFSFindPre) - length(findPos));
intTypeStore(2,2,fsfsTar) = intTypeStore(2,2,fsfsTar) + length(findPos);


%% Process correlation information: Look at FSFS -> FSFS inhibition before and after

%extract examples with no crossings before zero and more than 1 crossing
%after zero
findPre = find(alignedFSFSSig(:,3) >= 1);
findPost = find(alignedFSFSSig(:,2) >= 1);

combFSFSFindBoth = intersect(findPre,findPost);
%quality check, remove bad ones. 
qcCheck = zeros(length(combFSFSFindBoth),1);
for i = 1:length(combFSFSFindBoth)
    findNeg = find(bigSigNeg(alignedFSFS(combFSFSFindBoth(i)),:) < 0);
    tester = diff(findNeg);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

combFSFSFindBoth = combFSFSFindBoth(qcCheck == 1);
%determine number positive
findPos = find(alignedFSFSSig(combFSFSFindBoth,1) >1);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedFSFS(combFSFSFindBoth(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(4,1,fsfsTar) = intTypeStore(4,1,fsfsTar) + (length(combFSFSFindBoth) - length(findPos));
intTypeStore(4,2,fsfsTar) = intTypeStore(4,2,fsfsTar) + length(findPos);


%% Process correlation information: Look at FSFS -> FSFS no inhibition

findPre = find(alignedFSFSSig(:,3) <= 1);
findPost = find(alignedFSFSSig(:,2) <= 1);

combFSFSFindNone = intersect(findPre,findPost);

%determine number positive
findPos = find(alignedFSFSSig(combFSFSFindNone,1) >0);

%quality check, remove bad ones. 
qcCheck = zeros(length(findPos),1);
for i = 1:length(findPos)
    tarFind = find(bigSigPos(alignedFSFS(combFSFSFindNone(findPos(i))),:) > 0);
    tester = diff(tarFind);
    if min(tester) <= diffLim
        qcCheck(i) = 1;
    end
end

findPos = findPos(qcCheck == 1);

intTypeStore(1,1,fsfsTar) = intTypeStore(1,1,fsfsTar) + (length(combFSFSFindNone) - length(findPos));
intTypeStore(1,2,fsfsTar) = intTypeStore(1,2,fsfsTar) + length(findPos);

figure
for i = 1:14
subplot(3,5,i)
plot(bigxCorr(alignedFSFS(combFSFSFindBoth(i)),:))
end
















    
    
    
%% OLD STUFFs




%now, we need to find all examples where there is zero crossings before
%zero (L) and are crossings after zero (R)

    
    
    % I'm not longer certain this is a good way of approaching the issue.
    % This produces a bunch of disjointed datasets that are not necessarily
    % easy to put back together again. Instead, lets just march through
    % each array one by one and see what we pull out. We can log that
    % information then. 
    %lets see where we see significant correlation changes
    %start with shank 1, positive
    shank1UnitStorePos = zeros(1,2);
    shank1CorrStorePos = zeros(1,41);
    shank1CorrSigPos = zeros(1,41);
    shank1CorrShuffDatPos = zeros(3,41,1);
    counter = 1;
    for j = 1:size(corrData.SigCross1.PosWarn,1)
        %search along rows. row should be aligned to the target unit in
        %question.
        findSigPos = find(corrData.SigCross1.PosWarn(j,[1:j]) > 0);
        if findSigPos
            for k = 1:length(findSigPos)
                shank1UnitStorePos(counter,1) = j;
                shank1UnitStorePos(counter,2) = findSigPos(k);
                shank1CorrStorePos(counter,:) = squeeze(corrData.trueStoreShank1(j,findSigPos(k),:));
                shank1CorrSigPos = corrData.SigCross1.SigVals{j,findSigPos(k)};
                shank1CorrShuffDatPos(1,:,counter) = mean(corrData.ShuffStore1{j,findSigPos(k)}');
                shank1CorrShuffDatPos(2:3,:,counter) = prctile(corrData.ShuffStore1{j,findSigPos(k)},prctileBounds,2)';
                counter = counter + 1;
            end
        end
    end
    %now negative.
    
    %lets see where we see significant correlation changes
    %start with shank 1, positive
    shank1UnitStoreNeg = zeros(1,2);
    shank1CorrStoreNeg = zeros(1,41);
    shank1CorrSigNeg = zeros(1,41);
    shank1CorrShuffDatNeg = zeros(3,41,1);
    counter = 1;
    for j = 1:size(corrData.SigCross1.NegWarn,1)
        %search along rows. row should be aligned to the target unit in
        %question.
        findSigNeg = find(corrData.SigCross1.NegWarn(j,[1:j]) > 0);
        if findSigNeg
            for k = 1:length(findSigNeg)
                shank1UnitStoreNeg(counter,1) = j;
                shank1UnitStoreNeg(counter,2) = findSigNeg(k);
                shank1CorrStoreNeg(counter,:) = squeeze(corrData.trueStoreShank1(j,findSigNeg(k),:));
                shank1CorrSigNeg(counter,:) = corrData.SigCross1.SigVals{j,findSigNeg(k)};
                shank1CorrShuffDatNeg(1,:,counter) = mean(corrData.ShuffStore1{j,findSigNeg(k)}');
                shank1CorrShuffDatNeg(2:3,:,counter) = prctile(corrData.ShuffStore1{j,findSigNeg(k)},prctileBounds,2)';
                counter = counter + 1;
            end
        end
    end
    
    %now for shank 2. 
    shank2UnitStorePos = zeros(1,2);
    shank2CorrStorePos = zeros(1,41);
    shank2CorrSigPos = zeros(1,41);
    counter = 1;
    for j = 1:size(corrData.SigCross2.PosWarn,1)
        %search along rows. row should be aligned to the target unit in
        %question.
        findSigPos = find(corrData.SigCross2.PosWarn(j,[1:j]) > 0);
        if findSigPos
            for k = 1:length(findSigPos)
                shank2UnitStorePos(counter,1) = j;
                shank2UnitStorePos(counter,2) = findSigPos(k);
                shank2CorrStorePos(counter,:) = squeeze(corrData.trueStoreShank2(j,findSigPos(k),:));
                shank2CorrSigPos = corrData.SigCross2.SigVals{j,findSigPos(k)};
                counter = counter + 1;
            end
        end
    end
    %now negative.
    shank2UnitStoreNeg = zeros(1,2);
    shank2CorrStoreNeg = zeros(1,41);
    shank2CorrSigNeg = zeros(1,41);
    counter = 1;
    for j = 1:size(corrData.SigCross2.NegWarn,1)
        %search along rows. row should be aligned to the target unit in
        %question.
        findSigNeg = find(corrData.SigCross2.NegWarn(j,[1:j]) > 0);
        if findSigNeg
            for k = 1:length(findSigNeg)
                shank2UnitStoreNeg(counter,1) = j;
                shank2UnitStoreNeg(counter,2) = findSigNeg(k);
                shank2CorrStoreNeg(counter,:) = squeeze(corrData.trueStoreShank2(j,findSigNeg(k),:));
                shank2CorrSigNeg = corrData.SigCross2.SigVals{j,findSigNeg(k)};
                counter = counter + 1;
            end
        end
    end
    
    
    goodUnits = zeros(numUnits,1);
    %extract spike times, determine if exceeds RPV limits. 
    for j = 1:numUnits
        %go into spike count, get diff
        tmpISI = diff(s.(cluNames{j}).SpikeTimes);
        numSpikes = length(s.(cluNames{j}).SpikeTimes);
        tmpRPV = length(find(tmpISI < rpvCut))/numSpikes*100;
        rpvStore(j) = tmpRPV;
        if tmpRPV <= rpvLim
            goodUnits(j) = 1;
        else
            figure
            hist(tmpISI,rpvVect)
            xlim([0 0.02])
%             ylim([0 max(])
            title(num2str(tmpRPV))
        end
    end
    
    
end
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        