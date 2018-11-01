%This is code to do wrapper functions on analysisTuningWithWhite output
%data. Takes in s and masterData in, and is meant to output overall
%information about the population. 


%identify files, pull names, set up for loop for extraction

clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);


bigMaster = [];
bigMasterInd = 1;
respVect = [];
pkTroughRatioStore = [];
binValBigStore = [];
sigValBigStore = [];
latMapBigStore = [];
widthLatStore = [];
%actually extract files. 
for i = 1:numFiles
    %load the target file
    disp(strcat('Analyzing File:',num2str(i)))
    load(targetFiles{i})
    numUnits = size(masterData,1);
    numDBs = s.SoundData.NumDBs;
    pkTroughRatio = [];
    interpWaves = [];
    tempBinStore = [];
    tempSigStore = [];
    tempLatStore = [];
    tempMaxStore = [];
    tempHist = [];
    dbStore = [];
    freqStore = [];
    if numDBs == 5
        %pull from masterData, store in overall. 
        bigMaster(bigMasterInd:bigMasterInd + numUnits - 1,:) = masterData;

        %now pull overall designations of pos/neg/mix/no response
        [indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
        [indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');

        holder = masterData(:,[indPosSig,indNegSig]);
        holder(:,2) = holder(:,2) * -2;
        respVect(:,bigMasterInd:bigMasterInd + numUnits - 1) = (holder(:,1) + holder(:,2))'; %note that here, -2 = neg, -1 = mix, 0 = no, 1 = pos
        %find the peak to trough value ratio in waveforms. 
        numCells = length(s.DesignationName);
        desigName = s.DesignationName;
        for j = 1:numCells
            %pull up cell average waveforms
            cellWaves = s.(desigName{j}).AverageWaveForms;
            %now pull up which one is biggest
            waveMax = max(cellWaves);
            [maxVal maxInd] = max(waveMax);
            %now generate finely sampled version
            chosenWave = cellWaves(:,maxInd);
            interpVect = [1:0.1:40];
            interpWave = interp1(1:40,chosenWave,interpVect,'spline');
            %now find peak value

            [pkVal pkInd] = max(interpWave(100:end));
            pkInd = pkInd + 100 - 1;
            %now we need to find the minimum following the peak

            [troughVal troughInd] = min(interpWave(pkInd:end));
            troughInd = troughInd + pkInd - 1;
            pkTroughRatio(j) = pkVal/troughVal;
            interpWaves(j,:) = interpWave;

            %pull out binned values for entire tone period, as well as
            %significance
            tempBinStore(:,:,j) = s.(desigName{j}).BinTone;
            tempSigStore(:,:,j) = s.(desigName{j}).BinSigVals(:,:,2);
            tempLatStore(:,:,j) = s.(desigName{j}).LatencyMap;
            maxFind = find(s.(desigName{j}).PeakMapTone(2:end,:) == max(max(s.(desigName{j}).PeakMapTone(2:end,:))));
            if length(maxFind) == 1
                tempMaxStore(j) = maxFind;
            elseif length(maxFind) > 1
                maxFind = max(maxFind);
                tempMaxStore(j) = maxFind;
            else
                tempMaxStore(j) = 0;
            end
            %find appropriate histogram, store. 
            if length(maxFind) >= 1
                dbVal = floor((maxFind-1)/16);
                freqVal= mod(maxFind,16);
                tempHist(:,j) = squeeze(s.(desigName{j}).FreqDBHistograms(freqVal+1,dbVal+1,:));
                dbStore(j) = floor((maxFind-1)/16);
                freqStore(j) = mod(maxFind,16);
            else
                tempHist(:,j) = zeros(length(s.(desigName{j}).FreqDBHistograms(1,1,:)),1);
            end

        end
        pkTroughRatioStore(bigMasterInd:bigMasterInd + numUnits - 1) = pkTroughRatio;
        interpWaveStore(bigMasterInd:bigMasterInd + numUnits - 1,:) = interpWaves;
        binValBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempBinStore;
        sigValBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempSigStore;
        latMapBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempLatStore;
        widthLatStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = s.WidthLatData;
        bigMaxStore(bigMasterInd:bigMasterInd + numUnits - 1) = tempMaxStore;
        bigHistStore(:,bigMasterInd:bigMasterInd + numUnits - 1) = tempHist;
        bigDBStore(bigMasterInd:bigMasterInd + numUnits - 1) = dbStore;
        bigFreqStore(bigMasterInd:bigMasterInd + numUnits - 1) = freqStore;
        recStore(bigMasterInd:bigMasterInd + numUnits - 1) = i;
        
        %instead, lets just pull from posWidths and negWidths
        intFastPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.PosWidths((2:end),:,1));
        intFastNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NegWidths((2:end),:,1));
        intSlowPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.PosWidths((2:end),:,3));
        intSlowNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NegWidths((2:end),:,3));
        %store positive tuning widths
        widthStore(:,bigMasterInd:bigMasterInd + numUnits - 1,:) = s.PosWidths;
        %store BFs
        bfStore(bigMasterInd:bigMasterInd + numUnits - 1) = masterData(:,12);
        bigMasterInd = bigMasterInd + numUnits;
    elseif numDBs > 5
        %pull from masterData, store in overall. 
        bigMaster(bigMasterInd:bigMasterInd + numUnits - 1,:) = masterData;

        %now pull overall designations of pos/neg/mix/no response
        [indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
        [indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');

        holder = masterData(:,[indPosSig,indNegSig]);
        holder(:,2) = holder(:,2) * -2;
        respVect(:,bigMasterInd:bigMasterInd + numUnits - 1) = (holder(:,1) + holder(:,2))'; %note that here, -2 = neg, -1 = mix, 0 = no, 1 = pos
        %find the peak to trough value ratio in waveforms. 
        numCells = length(s.DesignationName);
        desigName = s.DesignationName;
        
        for j = 1:numCells
            %pull up cell average waveforms
            cellWaves = s.(desigName{j}).AverageWaveForms;
            %now pull up which one is biggest
            waveMax = max(cellWaves);
            [maxVal maxInd] = max(waveMax);
            %now generate finely sampled version
            chosenWave = cellWaves(:,maxInd);
            interpVect = [1:0.1:40];
            interpWave = interp1(1:40,chosenWave,interpVect,'spline');
            %now find peak value

            [pkVal pkInd] = max(interpWave(100:end));
            pkInd = pkInd + 100 - 1;
            %now we need to find the minimum following the peak

            [troughVal troughInd] = min(interpWave(pkInd:end));
            troughInd = troughInd + pkInd - 1;
            pkTroughRatio(j) = pkVal/troughVal;
            interpWaves(j,:) = interpWave;

            %pull out binned values for entire tone period, as well as
            %significance
            tempBinStore(:,:,j) = s.(desigName{j}).BinTone(:,end-5+1:end);
            tempSigStore(:,:,j) = s.(desigName{j}).BinSigVals(:,end-5+1:end,2);
            tempLatStore(:,:,j) = s.(desigName{j}).LatencyMap(:,end-5+1:end);
            maxFind = find(s.(desigName{j}).PeakMapTone(2:end,end-5+1:end) == max(max(s.(desigName{j}).PeakMapTone(2:end,end-5+1:end))));
            if length(maxFind) == 1
                tempMaxStore(j) = maxFind;
            elseif length(maxFind) > 1
                maxFind = max(maxFind);
                tempMaxStore(j) = max(maxFind);
            else
                tempMaxStore(j) = 0;
            end
            %find appropriate histogram, store. 
            if length(maxFind) >= 1
                dbVal = floor((maxFind-1)/16);
                freqVal= mod(maxFind,16);
                dbStore(j) = floor((maxFind-1)/16);
                freqStore(j) = mod(maxFind,16);
                tempHist(:,j) = squeeze(s.(desigName{j}).FreqDBHistograms(freqVal+1,dbVal+2,:));
            else
                tempHist(:,j) = zeros(length(s.(desigName{j}).FreqDBHistograms(1,1,:)),1);
            end
        end
        pkTroughRatioStore(bigMasterInd:bigMasterInd + numUnits - 1) = pkTroughRatio;
        interpWaveStore(bigMasterInd:bigMasterInd + numUnits - 1,:) = interpWaves;
        binValBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempBinStore;
        sigValBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempSigStore;
        latMapBigStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = tempLatStore;
        widthLatStore(:,:,bigMasterInd:bigMasterInd + numUnits - 1) = s.WidthLatData(:,end-5+1:end,:);
        bigMaxStore(bigMasterInd:bigMasterInd + numUnits - 1) = tempMaxStore;
        bigHistStore(:,bigMasterInd:bigMasterInd + numUnits - 1) = tempHist;
        bigDBStore(bigMasterInd:bigMasterInd + numUnits - 1) = dbStore;
        bigFreqStore(bigMasterInd:bigMasterInd + numUnits - 1) = freqStore;
        recStore(bigMasterInd:bigMasterInd + numUnits - 1) = i;
        %for "integral" of responses, we'll have to go into individual units,
        %and save both for short period and long period
    %     tempInd = 1;
    %     sigCutoff = 0.01;
    %     for j = 1:numUnits
    %         intFast(bigMasterInd + tempInd - 1) = length(find(squeeze(s.(s.DesignationName{j}).BinSigVals(2:end,:,1))<sigCutoff));
    %         intBig(bigMasterInd + tempInd - 1) = length(find(squeeze(s.(s.DesignationName{j}).BinSigVals(2:end,:,3))<sigCutoff));
    %         tempInd = tempInd + 1;
    %     end
        %instead, lets just pull from posWidths and negWidths
        intFastPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.PosWidths((2:end),:,1));
        intFastNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NegWidths((2:end),:,1));
        intSlowPos(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.PosWidths((2:end),:,3));
        intSlowNeg(:,bigMasterInd:bigMasterInd + numUnits - 1) = sum(s.NegWidths((2:end),:,3));
        %store positive tuning widths
        widthStore(:,bigMasterInd:bigMasterInd + numUnits - 1,:) = s.PosWidths(end-5+1:end,:,:);
        %store BFs
        bfStore(bigMasterInd:bigMasterInd + numUnits - 1) = masterData(:,12);
        bigMasterInd = bigMasterInd + numUnits;
    end
end



%now plot things!
[indCellType] = functionCellStringFind(masterHeader,'CellType');


%determine if there are units in each category
% findPVs = find(bigMaster(:,indCellType) == 1);
% 
% findMSNs = find(bigMaster(:,indCellType) == 0);

findCHATs = find(bigMaster(:,indCellType) == 2);


findPVs = find(bigMaster(:,5) < 0.0004);
findMSNs = find(bigMaster(:,5) > 0.0005);

%pull widths
pvWidths = widthStore(:,findPVs,:);
msnWidths = widthStore(:,findMSNs,:);
firstWidthMSN = msnWidths(:,:,2);
firstWidthPV = pvWidths(:,:,2);

%now we need to go through systematically. Store values like threshold
%value, width 10 dB above threshold, shape of response (monotonic
%increasing?)
infoStoreMSN = zeros(length(firstWidthMSN),6);
for i = 1:length(firstWidthMSN)
    %store the total number of DB steps with responses
    infoStoreMSN(i,1) = length(find(firstWidthMSN(:,i) > 0));
    %next, store the "integral" of responses, collapsing across dB
    infoStoreMSN(i,2) = sum(firstWidthMSN(:,i));
    %next, we need to try and look for threshold value. Do this by creeping
    %across. First, check to see if there are any values at all
    if length(find(firstWidthMSN(:,i) > 0)) <=1
        disp('Insufficient Responses Found')
        %report no threshold
        infoStoreMSN(i,3) = 0;
        %report no width
        infoStoreMSN(i,4) = 0;
        %report no shape
        infoStoreMSN(i,5) = 0;
    else
        disp('Sufficient Responses Found')
        %this assumes that there is stuff! now try and find things. 
        currVect = firstWidthMSN(:,i);
        %calculate threshold value in different ways, one by the first
        %point after, one at the last zero.
        threshVal1 = find(currVect > 0,1,'first');
        threshVal2 = find(currVect == 0,1,'last')+1;
        if threshVal1 == threshVal2 %if two thresholds are equal, then response is continuous
            threshVal = threshVal1;
            %now check to see if monotonic increasing
            infoStoreMSN(i,3) = threshVal; %store threshold value
            %now store width at level above threshold. If threshold equals
            %lowest value, take that width. 
            if threshVal < length(currVect) - 1;
                infoStoreMSN(i,4) = currVect(threshVal+1);
            else
                infoStoreMSN(i,4) = currVect(threshVal);
            end
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStoreMSN(i,5) = 1; %store monotonic increasing shape
            else
                infoStoreMSN(i,5) = 2; %store value as non-monotonic
            end
        elseif threshVal1 == 1 & isempty(threshVal2)
            infoStoreMSN(i,3) = 1;%store threshold as 1
            infoStoreMSN(i,4) = currVect(2); %store width
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStoreMSN(i,5) = 1; %store monotonic increasing shape
            else
                infoStoreMSN(i,5) = 2; %store value as non-monotonic
            end
        else
            threshVal = threshVal2;
            infoStoreMSN(i,3) = threshVal; %store threshold value. This is still useful information
            infoStoreMSN(i,4) = 0;
            infoStoreMSN(i,5) = 2;
        end
    end
end

%do the same for PVs
infoStorePV = zeros(length(firstWidthPV),6);
for i = 1:length(firstWidthPV)
    %store the total number of DB steps with responses
    infoStorePV(i,1) = length(find(firstWidthPV(:,i) > 0));
    %next, store the "integral" of responses, collapsing across dB
    infoStorePV(i,2) = sum(firstWidthPV(:,i));
    %next, we need to try and look for threshold value. Do this by creeping
    %across. First, check to see if there are any values at all
    if length(find(firstWidthPV(:,i) > 0)) <=1
        disp('Insufficient Responses Found')
        %report no threshold
        infoStorePV(i,3) = 0;
        %report no width
        infoStorePV(i,4) = 0;
        %report no shape
        infoStorePV(i,5) = 0;
    else
        disp('Sufficient Responses Found')
        %this assumes that there is stuff! now try and find things. 
        currVect = firstWidthPV(:,i);
        %calculate threshold value in different ways, one by the first
        %point after, one at the last zero.
        threshVal1 = find(currVect > 0,1,'first');
        threshVal2 = find(currVect == 0,1,'last')+1;
        if threshVal1 == threshVal2 %if two thresholds are equal, then response is continuous
            threshVal = threshVal1;
            %now check to see if monotonic increasing
            infoStorePV(i,3) = threshVal; %store threshold value
            if threshVal < length(currVect) - 1;
                infoStorePV(i,4) = currVect(threshVal+1);
            else
                infoStorePV(i,4) = currVect(threshVal);
            end
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStorePV(i,5) = 1; %store monotonic increasing shape
            else
                infoStorePV(i,5) = 2; %store value as non-monotonic
            end
        elseif threshVal1 == 1 & isempty(threshVal2)
            infoStorePV(i,3) = 1;%store threshold as 1
            infoStorePV(i,4) = currVect(2); %store width
            
            vectDiff = diff(currVect);
            findNeg = find(vectDiff < 0);
            if isempty(findNeg) %no negative changes! yay!
                infoStorePV(i,5) = 1; %store monotonic increasing shape
            else
                infoStorePV(i,5) = 2; %store value as non-monotonic
            end
        else
            threshVal = threshVal2;
            infoStorePV(i,3) = threshVal; %store threshold value. This is still useful information
            infoStorePV(i,4) = 0;
            infoStorePV(i,5) = 2;
        end
    end
end

condensedMSN = infoStoreMSN(infoStoreMSN(:,5) == 1,:);

condensedPV = infoStorePV(infoStorePV(:,5) == 1,:);



[indPkTr] = functionCellStringFind(masterHeader,'PeakTrough');
[indISI] = functionCellStringFind(masterHeader,'isiCov');
[indBaseFire] = functionCellStringFind(masterHeader,'BaselineRate');
[indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
[indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');

%make vector for width. if neuron responds to everything, that makes 5 *
%16, or 80. 
widthHistVect = [0:1:40];

%% NOW LETS DO PLOTTING

%% first figure, looking at broad general things, like numbers of cells and firing rates
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])

%column 1
%plot spike width vs coefficient of variation
subplot(4,4,1)
hold on

plot(bigMaster(bigMaster(:,indCellType) == 0,indPkTr),bigMaster(bigMaster(:,indCellType) == 0,indISI),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),bigMaster(bigMaster(:,indCellType) == 1,indISI),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),bigMaster(bigMaster(:,indCellType) == 2,indISI),'g.')
xlabel('Peak Trough (ms)')
ylabel('ISI Coefficient of Variation')

%now plot out spike width vs peak trough ratio
subplot(4,4,5)
hold on

plot(bigMaster(:,indPkTr),pkTroughRatioStore(:),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),pkTroughRatioStore(bigMaster(:,indCellType) == 1),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),pkTroughRatioStore(bigMaster(:,indCellType) == 2),'g.')
xlabel('Peak Trough (ms)')
ylabel('Peak Trough Amplitude Ratio')

%plot out spike width with FR
subplot(4,4,9)
hold on

plot(bigMaster(:,indPkTr),bigMaster(:,indBaseFire),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),bigMaster(bigMaster(:,indCellType) == 1,indBaseFire),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),bigMaster(bigMaster(:,indCellType) == 2,indBaseFire),'g.')
xlabel('Peak Trough (ms)')
ylabel('Baseline Firing Rate')

%plot out pie chart of difference cells
subplot(4,4,13)
cellDist = [length(findMSNs),length(findPVs),length(findCHATs)];
pie(cellDist)
labels = {strcat('MSNs(',num2str(length(findMSNs)),')'),strcat('PVs(',num2str(length(findPVs)),')'),strcat('ChATs(',num2str(length(findCHATs)),')')};
detZero = find(cellDist == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')

%column 2, plot out pie charts of MSNs

subplot(4,4,2)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findMSNs,:);
plot(waveHolder')
title('MSN Waveforms')

%plot distribution of response types
subplot(4,4,6)
holder = bigMaster(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('MSNs n=',num2str(length(findMSNs))))

%plot out integrals of response plots
subplot(4,4,10)
holder = intFastPos(findMSNs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))


subplot(4,4,14)
hold on
holder = intSlowPos(findMSNs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))


%column 3, plot out PV stuff

subplot(4,4,3)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findPVs,:);
plot(waveHolder')
title('PV Waveforms')

%plot distribution of response types
subplot(4,4,7)
holder = bigMaster(findPVs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('PVs n=',num2str(length(findPVs))))

%plot out integrals of response plots
subplot(4,4,11)
holder = intFastPos(findPVs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

subplot(4,4,15)
hold on
holder = intSlowPos(findPVs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

%column 4, plot out pie charts of CHATs

subplot(4,4,4)
hold on
waveHolder = [];
waveHolder = interpWaveStore(findCHATs,:);
plot(waveHolder')
title('CHAT Waveforms')

%plot distribution of response types
subplot(4,4,8)
holder = bigMaster(findCHATs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','eastoutside','Orientation','vertical')
title(strcat('CHATs n=',num2str(length(findCHATs))))

%plot out integrals of response plots
subplot(4,4,12)
holder = intFastPos(findCHATs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

subplot(4,4,16)
hold on
holder = intSlowPos(findCHATs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))


spikeGraphName = 'WrapperFigure1Overview';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')



%lets try plotting out binned responses
% hFig = figure;
% set(hFig, 'Position', [10 80 1900 1000])
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% for i = 1:200
%     subplot(10,20,i)
%     imagesc(binValBigStore(:,:,i)')
%     colormap('parula')
%     set(gca,'YTick',[]);
%     set(gca,'XTick',[]);
%     set(gca,'Ydir','reverse')
% end
% 
% %plot without white noise responses
% hFig = figure;
% set(hFig, 'Position', [10 80 1900 1000])
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% for i = 1:length(findMSNs)
%     subplot(15,20,i)
%     imagesc(binValBigStore(2:end,:,i)')
%     colormap('parula')
%     set(gca,'YTick',[]);
%     set(gca,'XTick',[]);
%     set(gca,'Ydir','reverse')
% end


%how lets look at sigValBigStore
sigValConv = sigValBigStore;
sigValConv(sigValConv <= 0.001) = 4;
sigValConv(sigValConv <= 0.01) = 3;
sigValConv(sigValConv <= 0.05) = 2;
sigValConv(sigValConv <= 1) = 1;
% 
% hFig = figure;
% set(hFig, 'Position', [10 80 1900 1000])
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
% for i = 1:length(findMSNs)
%     subplot(15,20,i)
%     imagesc(sigValConv(2:end,:,i)')
%     colormap('parula')
%     set(gca,'YTick',[]);
%     set(gca,'XTick',[]);
%     set(gca,'Ydir','reverse')
% end

clims = [1 3];

%% Plot out tuning curves without white noise
%just do PVs
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findPVs)
    subplot(8,6,i)
    imagesc(binValBigStore(2:end,:,findPVs(i))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'pvBinStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findPVs)
    subplot(8,6,i)
    imagesc(sigValConv(2:end,:,findPVs(i))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end
spikeGraphName = 'pvSigStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%and just MSNs
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findMSNs)
    subplot(15,20,i)
    imagesc(binValBigStore(2:end,:,findMSNs(i))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'msnBinStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findMSNs)
    subplot(15,20,i)
    imagesc(sigValConv(2:end,:,findMSNs(i))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'msnSigStoreTuning';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

% based on this, it seems like there is something of a difference between
% FSIs and MSNs

%% lets try plotting just white noise responses
%just do PVs

clims = [1 3];

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findPVs)
    subplot(1,length(findPVs),i)
    imagesc(binValBigStore(1,:,findPVs(i))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findPVs)
    subplot(1,length(findPVs),i)
    imagesc(sigValConv(1,:,findPVs(i))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

%and just MSNs
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findMSNs)
    subplot(15,20,i)
    imagesc(binValBigStore(1,:,findMSNs(i))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(findMSNs)
    subplot(15,20,i)
    imagesc(sigValConv(1,:,findMSNs(i))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end


%% now lets look at latency
%180905 new attempt. Will threshold by units that respond to at least two
%things during tone period significantly. 
tester = max(widthStore(:,:,2));%###########
tarCells = find(tester >= 3);

%first, we want to eliminate non-significant latencies
sigMin = 0.01; %minimum significance value
binSig = sigValBigStore;
binSig(binSig > sigMin) = NaN;
binSig(binSig <= sigMin) = 1;

latConv = binSig.*latMapBigStore;
widthLatConv = binSig.*widthLatStore;
%now delete zero latency values
latConv(latConv == 0) = NaN;

%lets separate white noise from the rest
latConvWhite = latConv(1,:,:);
latConvTone = latConv(2:end,:,:);
latConvWidthTone = widthLatConv(2:end,:,:);

%now lets extract targeted latencies. 
for i = 1:length(tarCells);
    tempLat = latConvTone(:,:,tarCells(i));
    tarLats(i) = tempLat(bigMaxStore(tarCells(i)));
end
%QC check, remove all latencies and tarCells that are NaNs
nanFind = isnan(tarLats);
tarLats(nanFind) = [];
tarCells(nanFind) = [];

[C,tarPVs,ib] = intersect(tarCells,findPVs);
[C,tarMSNs,ib] = intersect(tarCells,findMSNs);

latHistVect = [0:0.001:0.1];

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(3,2,1)
hist(tarLats(tarPVs),latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(tarPVs)), ' FSI Min Latency Tone'))

% subplot(3,2,2)
% hist(minLatTonePV,latHistVect)
% xlim([latHistVect(1) latHistVect(end)])
% title(strcat(num2str(length(find(~isnan(minLatTonePV)))),' FSI Min Latency Pure Tone'))
% % title('FSI Min Latency Pure Tone')

subplot(3,2,3)
hist(tarLats(tarMSNs),latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(tarMSNs)),' MSN Min Latency Tone'))
% title('MSN Min Latency White Noise')

% subplot(3,2,4)
% hist(minLatToneMSN,latHistVect)
% xlim([latHistVect(1) latHistVect(end)])
% title(strcat(num2str(length(find(~isnan(minLatToneMSN)))),' MSN Min Latency Pure Tone'))
% % title('MSN Min Latency Pure Tone')

subplot(3,1,3)
hold on
bar(1:2,[nanmean(tarLats(tarPVs)),nanmean(tarLats(tarMSNs))],'w')
errorbar(1:2,[nanmean(tarLats(tarPVs)),nanmean(tarLats(tarMSNs))],[nanstd(tarLats(tarPVs)),nanstd(tarLats(tarMSNs))])
% xticks([1:4])
% xticklabels({'FSI White','MSN White','FSI Tone','MSN Tone'})

spikeGraphName = 'latencyPlotSelected';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%the results of this demonstrate that there is a weak difference in latency
%that is insignificant. Looking at histograms, looks like this might be in
%part due to a number of slower responding PV cells. 

%lets try and look on a recording by recording basis
tarRec = recStore(tarCells)';
tarRec(:,2) = 0;
tarRec(tarPVs,2) = 1;
tarRec(:,3) = 0;
%check to see if specific recording has both PVs and MSNs
for i = 1:numFiles
    tempStore = tarRec(tarRec(:,1) == i,2);
    if sum(tempStore) == 0 | sum(tempStore) == length(tempStore)
        disp('Only One Cell Type')
    else
        disp('Multiple Cell Types')
        tarRec(tarRec(:,1) == i,3) = i;
    end
end

tarRec(:,4) = tarLats;

figure
hold on
plot(mean(bigHistStore(:,tarCells(tarPVs))'))
plot(mean(bigHistStore(:,tarCells(tarMSNs))'),'r')

%Doesnt look like a per-recording analysis will pull out anything different
%really. 


minLatWhitePV = squeeze(min(latConvWhite(:,:,findPVs)));
minLatWhiteMSN = squeeze(min(latConvWhite(:,:,findMSNs)));
minLatTonePV = squeeze(min(min(latConvTone(:,:,findPVs))));
minLatToneMSN = squeeze(min(min(latConvTone(:,:,findMSNs))));

minLatWhitePVWidth = squeeze(min(latConvWidthTone(1,:,findPVs)));
minLatWhiteMSNWidth = squeeze(min(latConvWidthTone(1,:,findMSNs)));
minLatTonePVWidth = squeeze(min(min(latConvWidthTone(:,:,findPVs))));
minLatToneMSNWidth = squeeze(min(min(latConvWidthTone(:,:,findMSNs))));

latHistVect = [0:0.001:0.1];

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(3,2,1)
hist(minLatWhitePV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatWhitePV)))), ' FSI Min Latency White Noise'))

subplot(3,2,2)
hist(minLatTonePV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatTonePV)))),' FSI Min Latency Pure Tone'))
% title('FSI Min Latency Pure Tone')

subplot(3,2,3)
hist(minLatWhiteMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatWhiteMSN)))),' MSN Min Latency White Noise'))
% title('MSN Min Latency White Noise')

subplot(3,2,4)
hist(minLatToneMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatToneMSN)))),' MSN Min Latency Pure Tone'))
% title('MSN Min Latency Pure Tone')

subplot(3,1,3)
hold on
bar(1:4,[nanmean(minLatWhitePV),nanmean(minLatWhiteMSN),nanmean(minLatTonePV),nanmean(minLatToneMSN)],'w')
errorbar(1:4,[nanmean(minLatWhitePV),nanmean(minLatWhiteMSN),nanmean(minLatTonePV),nanmean(minLatToneMSN)],[nanstd(minLatWhitePV),nanstd(minLatWhiteMSN),nanstd(minLatTonePV),nanstd(minLatToneMSN)])
% xticks([1:4])
% xticklabels({'FSI White','MSN White','FSI Tone','MSN Tone'})

spikeGraphName = 'latencyPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%plot using latencies from width finder. 
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(3,2,1)
hist(minLatWhitePVWidth,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatWhitePVWidth)))), ' FSI Min Latency White Noise'))

subplot(3,2,2)
hist(minLatTonePVWidth,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatTonePVWidth)))),' FSI Min Latency Pure Tone'))
% title('FSI Min Latency Pure Tone')

subplot(3,2,3)
hist(minLatWhiteMSNWidth,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatWhiteMSNWidth)))),' MSN Min Latency White Noise'))
% title('MSN Min Latency White Noise')

subplot(3,2,4)
hist(minLatToneMSNWidth,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatToneMSNWidth)))),' MSN Min Latency Pure Tone'))
% title('MSN Min Latency Pure Tone')

subplot(3,1,3)
hold on
bar(1:4,[nanmean(minLatWhitePVWidth),nanmean(minLatWhiteMSNWidth),nanmean(minLatTonePVWidth),nanmean(minLatToneMSNWidth)],'w')
errorbar(1:4,[nanmean(minLatWhitePVWidth),nanmean(minLatWhiteMSNWidth),nanmean(minLatTonePVWidth),nanmean(minLatToneMSNWidth)],[nanstd(minLatWhitePVWidth),nanstd(minLatWhiteMSNWidth),nanstd(minLatTonePVWidth),nanstd(minLatToneMSNWidth)])
% xticks([1:4])
% xticklabels({'FSI White','MSN White','FSI Tone','MSN Tone'})

spikeGraphName = 'WidLatencyPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')
%based on this, there clearly is some failure of width based latency
%calculations. 


%now lets look at specific amplitudes?

minLatWhitePV = squeeze((latConvWhite(:,5,findPVs)));
minLatWhiteMSN = squeeze((latConvWhite(:,5,findMSNs)));
minLatTonePV = squeeze((min(latConvTone(:,5,findPVs))));
minLatToneMSN = squeeze((min(latConvTone(:,5,findMSNs))));

latHistVect = [0:0.001:0.1];

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(3,2,1)
hist(minLatWhitePV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatWhitePV)))), ' FSI Min Latency White Noise'))

subplot(3,2,2)
hist(minLatTonePV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatTonePV)))),' FSI Min Latency Pure Tone'))
% title('FSI Min Latency Pure Tone')

subplot(3,2,3)
hist(minLatWhiteMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatWhiteMSN)))),' MSN Min Latency White Noise'))
% title('MSN Min Latency White Noise')

subplot(3,2,4)
hist(minLatToneMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minLatToneMSN)))),' MSN Min Latency Pure Tone'))
% title('MSN Min Latency Pure Tone')

subplot(3,1,3)
hold on
bar(1:4,[nanmean(minLatWhitePV),nanmean(minLatWhiteMSN),nanmean(minLatTonePV),nanmean(minLatToneMSN)],'w')
errorbar(1:4,[nanmean(minLatWhitePV),nanmean(minLatWhiteMSN),nanmean(minLatTonePV),nanmean(minLatToneMSN)],[nanstd(minLatWhitePV),nanstd(minLatWhiteMSN),nanstd(minLatTonePV),nanstd(minLatToneMSN)])
% xticks([1:4])
% xticklabels({'FSI White','MSN White','FSI Tone','MSN Tone'})

spikeGraphName = 'latencyPlot70DB';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')




%% lets try plotting out only the "classic" responses 

pvTarget = find(infoStorePV(:,5)==1);
msnTarget = find(infoStoreMSN(:,5)==1);
%limit by same 5 value limit. 

sigLim = 5;
%find 
tester = binSig(2:end,:,findPVs);
counter = 1;
for i = 1:length(findPVs)
    sumVal = sum(nansum(tester(:,:,i)));
    if sumVal >= sigLim
        sigPVs(counter) = i;
        counter = counter + 1;
    end
end

tester = binSig(2:end,:,findMSNs);
counter = 1;
for i = 1:length(findMSNs)
    sumVal = sum(nansum(tester(:,:,i)));
    if sumVal >= sigLim
        sigMSNs(counter) = i;
        counter = counter + 1;
    end
end

pvTarget = intersect(pvTarget,sigPVs);
msnTarget = intersect(msnTarget,sigMSNs);


clims = [1 3];
%just do PVs
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(pvTarget)
    subplot(5,3,i)
    imagesc(binValBigStore(2:end,:,findPVs(pvTarget(i)))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'pvBinStoreTuningSelected';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(pvTarget)
    subplot(5,3,i)
    imagesc(sigValConv(2:end,:,findPVs(pvTarget(i)))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end
spikeGraphName = 'pvSigStoreTuningSelected';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%and just MSNs
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(msnTarget)
    subplot(6,5,i)
    imagesc(binValBigStore(2:end,:,findMSNs(msnTarget(i)))')
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'msnBinStoreTuningSelected';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.005 0.005], [0.005 0.005]);
for i = 1:length(msnTarget)
    subplot(6,5,i)
    imagesc(sigValConv(2:end,:,findMSNs(msnTarget(i)))',clims)
    colormap('parula')
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca,'Ydir','reverse')
end

spikeGraphName = 'msnSigStoreTuningSelected';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')




%% Now lets look at what the threshold of activation is. 

%first do PVs. 
threshStorePV = [];
for i = 1:length(pvTarget)
    tester = sigValConv(2:end,:,findPVs(pvTarget(i)));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    threshStorePV(i) = testFind;
end

%now do MSNs
threshStoreMSN = [];
for i = 1:length(msnTarget)
    tester = sigValConv(2:end,:,findMSNs(msnTarget(i)));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    threshStoreMSN(i) = testFind;
end


%plot this out. 
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,2,1)
hist(threshStorePV,[1:1:5])
xlim([0.5 5.5])
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});
title(strcat('Amplitude Threshold (FSIs) Mean:',num2str(mean(threshStorePV))))

subplot(2,2,2)
hist(threshStoreMSN,[1:1:5])
xlim([0.5 5.5])
title(strcat('Amplitude Threshold (MSNs) Mean:',num2str(mean(threshStoreMSN))))
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});

subplot(2,2,3)
hist(firstWidthPV(5,pvTarget),[1:1:16])
xlim([0.5 16.5])
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('Tuning Width (FSIs) Mean:',num2str(mean(firstWidthPV(5,pvTarget)))))
% title('Tuning Width (FSIs)')

subplot(2,2,4)
hist(firstWidthMSN(5,msnTarget),[1:1:16])
xlim([0.5 16.5])
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('Tuning Width (MSNs) Mean:',num2str(mean(firstWidthMSN(5,msnTarget)))))

spikeGraphName = 'SelectedIntensityThresholdAndWidthPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

ranksum(threshStorePV,threshStoreMSN)
ranksum(firstWidthPV(5,pvTarget),firstWidthMSN(5,msnTarget))

%plot for everything. 
%first do PVs. 
threshStorePV = [];
for i = 1:length(findPVs)
    tester = sigValConv(2:end,:,findPVs(i));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    if isempty(testFind)
        threshStorePV(i) = NaN;
    else
        threshStorePV(i) = testFind;
    end
    
end

%now do MSNs
threshStoreMSN = [];
for i = 1:length(findMSNs)
    tester = sigValConv(2:end,:,findMSNs((i)));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    if isempty(testFind)
        threshStoreMSN(i) = NaN;
    else
        threshStoreMSN(i) = testFind;
    end
%     threshStoreMSN(i) = testFind;
end

%plot this out. 
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,2,1)
hist(threshStorePV,[1:1:7])
xlim([0.5 5.5])
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});
title(strcat('Amplitude Threshold (FSIs) Mean:',num2str(nanmean(threshStorePV))))

subplot(2,2,2)
hist(threshStoreMSN,[1:1:7])
xlim([0.5 5.5])
title(strcat('Amplitude Threshold (MSNs) Mean:',num2str(nanmean(threshStoreMSN))))
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});

subplot(2,2,3)
hist(firstWidthPV(5,:),[1:1:16])
xlim([0.5 16.5])
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('Tuning Width (FSIs) Mean:',num2str(mean(firstWidthPV(5,:)))))
% title('Tuning Width (FSIs)')

subplot(2,2,4)
hist(firstWidthMSN(5,:),[1:1:16])
xlim([0.5 16.5])
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('Tuning Width (MSNs) Mean:',num2str(mean(firstWidthMSN(5,:)))))

spikeGraphName = 'GeneralIntensityThresholdAndWidthPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

ranksum(threshStorePV,threshStoreMSN)
ranksum(firstWidthPV(5,pvTarget),firstWidthMSN(5,msnTarget))

%These plots and stats indicate the following: 
%overall, if we look at all cells, there isnt a significant change in
%threshold between MSNs and FSIs. This is code the relies upon the first
%detected significant response, however, so it is susceptible to noise.
%However, things do look somewhat different in the plot.
%Looking at the overall though, there definitely appears to be a
%significant difference in the tuning width, which is in part driven
%strongly by the fact that many MSNs are non-responsive. 

%If we look at only neurons that have somewhat correctly shaped tuning
%curves, we find that both threshold and width are skewed in favor of FSIs
%in a significant fashion. 

%as a next step, lets examine the plots with MSNs taht are labeled as
%having positive responses. 

%plot for JUST significant responders
%first do PVs. 

sigLim = 5;
%find 
tester = binSig(2:end,:,findPVs);
counter = 1;
for i = 1:length(findPVs)
    sumVal = sum(nansum(tester(:,:,i)));
    if sumVal >= sigLim
        sigPVs(counter) = i;
        counter = counter + 1;
    end
end

tester = binSig(2:end,:,findMSNs);
counter = 1;
for i = 1:length(findMSNs)
    sumVal = sum(nansum(tester(:,:,i)));
    if sumVal >= sigLim
        sigMSNs(counter) = i;
        counter = counter + 1;
    end
end


threshStorePV = [];
for i = 1:length(sigPVs)
    tester = sigValConv(2:end,:,findPVs(sigPVs(i)));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    if isempty(testFind)
        threshStorePV(i) = NaN;
    else
        threshStorePV(i) = testFind;
    end
    
end

%now do MSNs
threshStoreMSN = [];
for i = 1:length(sigMSNs)
    tester = sigValConv(2:end,:,findMSNs(sigMSNs(i)));
    condTester = max(tester);
    %find first value == 3
    testFind = find(condTester >= 3,1,'first');
    if isempty(testFind)
        threshStoreMSN(i) = NaN;
    else
        threshStoreMSN(i) = testFind;
    end
%     threshStoreMSN(i) = testFind;
end

%plot this out. 
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.07 0.07], [0.07 0.07]);
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,2,1)
hist(threshStorePV,[1:1:7])
xlim([0.5 5.5])
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});
title(strcat('Amplitude Threshold (FSIs) Mean:',num2str(nanmean(threshStorePV))))

subplot(2,2,2)
hist(threshStoreMSN,[1:1:7])
xlim([0.5 5.5])
title(strcat('Amplitude Threshold (MSNs) Mean:',num2str(nanmean(threshStoreMSN))))
ylabel('Number of Cells')
xlabel('Amplitude of Threshold Response (dB)')
set(gca,'XTickLabel',{'<=20','30','40','60','70'});

subplot(2,2,3)
hist(firstWidthPV(5,(sigPVs)),[1:1:16])
xlim([0.5 16.5])
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('Tuning Width (FSIs) Mean:',num2str(mean(firstWidthPV(5,sigPVs)))))
% title('Tuning Width (FSIs)')

subplot(2,2,4)
hist(firstWidthMSN(5,(sigMSNs)),[1:1:16])
xlim([0.5 16.5])
ylabel('Number of Cells')
xlabel('Tuning Curve Width (0.2 octaves)')
title(strcat('Tuning Width (MSNs) Mean:',num2str(mean(firstWidthMSN(5,sigMSNs)))))

spikeGraphName = 'val5IntensityThresholdAndWidthPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

ranksum(threshStorePV,threshStoreMSN)
ranksum(firstWidthPV(5,sigPVs),firstWidthMSN(5,sigMSNs))

%If we only include neurons that have responses that are not to white
%noises at at least 5 points, then we see a significant difference in
%threshold and width. 

%lets apply these groups to latency calculations

minLatVals = squeeze(min(latConvTone));

%lets get values for 5 response groups
minLat5ValPV = minLatVals(:,findPVs(sigPVs));
minLat5ValMSN = minLatVals(:,findMSNs(sigMSNs));
%get values for "good looking" units
minLatTarPV = minLatVals(:,findPVs(pvTarget));
minLatTarMSN = minLatVals(:,findMSNs(msnTarget));

minmin5ValPV = min(minLat5ValPV);
minmin5ValMSN = min(minLat5ValMSN);

minminLatTarPV = min(minLatTarPV);
minminLatTarMSN = min(minLatTarMSN);


hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(2,2,1)
hist(minmin5ValPV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minmin5ValPV)))), ' FSI 5Val Min Lat All DB'))

subplot(2,2,2)
hist(minminLatTarPV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minminLatTarPV)))), ' FSI Tar Min Lat All DB'))

subplot(2,2,3)
hist(minmin5ValMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minmin5ValMSN)))),' MSN Min Lat All DB'))
% title('FSI Min Latency Pure Tone')

subplot(2,2,4)
hist(minminLatTarMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minminLatTarMSN)))),' MSN Min Lat All DB'))
% title('FSI Min Latency Pure Tone')

spikeGraphName = 'AllDB5Val&TarPVandMSNLatencies';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

nanmean(minmin5ValPV)
nanmean(minmin5ValMSN)
nanmean(minminLatTarPV)
nanmean(minminLatTarMSN)

ranksum(minmin5ValPV,minmin5ValMSN)
ranksum(minminLatTarPV,minminLatTarMSN)


%now lets try just looking at max amplitude. 

minmin5ValPV = (minLat5ValPV(5,:));
minmin5ValMSN = (minLat5ValMSN(5,:));

minminLatTarPV = (minLatTarPV(5,:));
minminLatTarMSN = (minLatTarMSN(5,:));

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(2,2,1)
hist(minmin5ValPV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minmin5ValPV)))), ' FSI 5Val Min Lat 70 DB'))

subplot(2,2,2)
hist(minminLatTarPV,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minminLatTarPV)))), ' FSI Tar Min Lat 70 DB'))

subplot(2,2,3)
hist(minmin5ValMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minmin5ValMSN)))),' MSN Min Lat 70 DB'))
% title('FSI Min Latency Pure Tone')

subplot(2,2,4)
hist(minminLatTarMSN,latHistVect)
xlim([latHistVect(1) latHistVect(end)])
title(strcat(num2str(length(find(~isnan(minminLatTarMSN)))),' MSN Min Lat 70 DB'))
% title('FSI Min Latency Pure Tone')

spikeGraphName = '70DB5Val&TarPVandMSNLatencies';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')




nanmean(minmin5ValPV)
nanmean(minmin5ValMSN)
nanmean(minminLatTarPV)
nanmean(minminLatTarMSN)

ranksum(minmin5ValPV,minmin5ValMSN)
ranksum(minminLatTarPV,minminLatTarMSN)




%% now lets look at best frequencies!

bfs = bigMaster(:,12);
BFpv = bfs(findPVs);
BFmsn = bfs(findMSNs);

calibChart = [4000	2	0.7943282347
4287.09385	0	1
4594.79342	0.5	0.9440608763
4924.577653	1	0.8912509381
5278.031643	1.5	0.8413951416
5656.854249	1.5	0.8413951416
6062.866266	2.2	0.7762471166
6498.019171	1.5	0.8413951416
6964.404506	2	0.7943282347
7464.263932	4.5	0.5956621435
8000	6.2	0.4897788194
8574.1877	4.3	0.6095368972
9189.58684	3.6	0.660693448
9849.155307	6.3	0.4841723676
10556.06329	5.1	0.5559042573
11313.7085	3.8	0.645654229
12125.73253	1.5	0.8413951416
12996.03834	3.5	0.6683439176
13928.80901	3.7	0.6531305526
14928.52786	3.8	0.645654229
16000	3.5	0.6683439176
17148.3754	2.5	0.7498942093
18379.17368	2	0.7943282347
19698.31061	6.2	0.4897788194
21112.12657	7.8	0.4073802778
22627.417	8.75	0.3651741273
24251.46506	10	0.316227766
25992.07668	14	0.1995262315
27857.61803	12	0.2511886432
29857.05573	13.2	0.2187761624
32000	15.6	0.1659586907
34296.7508	16.5	0.1496235656
36758.34736	18.2	0.1230268771
39396.62123	20	0.1
42224.25314	16	0.1584893192
45254.834	18.7	0.1161448614
48502.93013	14.2	0.19498446
51984.15337	15.5	0.1678804018
55715.23605	11.7	0.2600159563
59714.11146	11.7	0.2600159563
64000	10	0.316227766];

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
subplot(1,2,1)
hist(BFpv,calibChart(:,1))
set(gca,'xscale','log') 
xlim([4000 32000])
title('Best Frequencies PV')

subplot(1,2,2)
hist(BFmsn,calibChart(:,1))
set(gca,'xscale','log') 
xlim([4000 32000])
title('Best Frequencies MSN')


spikeGraphName = 'histOfBFsPVMSN';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')












