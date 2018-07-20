%This is code to do wrapper functions on analysisTuningWithWhite output
%data. Takes in s and masterData in, and is meant to output overall
%information about the population. 


%identify files, pull names, set up for loop for extraction

% clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);


bigMaster = [];
bigMasterInd = 1;
respVect = [];

%actually extract files. 
for i = 1:numFiles
    %load the target file
    load(targetFiles{i})
    numUnits = size(masterData,1);
    
    %pull from masterData, store in overall. 
    bigMaster(bigMasterInd:bigMasterInd + numUnits - 1,:) = masterData;
    
    %now pull overall designations of pos/neg/mix/no response
    [indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
    [indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');
    
    holder = masterData(:,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    respVect(:,bigMasterInd:bigMasterInd + numUnits - 1) = (holder(:,1) + holder(:,2))'; %note that here, -2 = neg, -1 = mix, 0 = no, 1 = pos
    
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
    widthStore(:,bigMasterInd:bigMasterInd + numUnits - 1,:) = s.PosWidths;
    %store BFs
    bfStore(bigMasterInd:bigMasterInd + numUnits - 1) = masterData(:,12);
    bigMasterInd = bigMasterInd + numUnits;
end



%now plot things!
[indCellType] = functionCellStringFind(masterHeader,'CellType');


%determine if there are units in each category
findPVs = find(bigMaster(:,indCellType) == 1);

findMSNs = find(bigMaster(:,indCellType) == 0);

findCHATs = find(bigMaster(:,indCellType) == 2);

%pull widths
pvWidths = widthStore(:,findPVs,:);
msnWidths = widthStore(:,findMSNs,:);
firstWidthMSN = msnWidths(:,:,1);
firstWidthPV = pvWidths(:,:,1);

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
[indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
[indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');

%make vector for width. if neuron responds to everything, that makes 5 *
%16, or 80. 
widthHistVect = [0:1:80];


subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
%% Column 1
%plot spike width vs coefficient of variation
subplot(3,4,1)
hold on

plot(bigMaster(:,indPkTr),bigMaster(:,indISI),'k.')
plot(bigMaster(bigMaster(:,indCellType) == 1,indPkTr),bigMaster(bigMaster(:,indCellType) == 1,indISI),'r.')
plot(bigMaster(bigMaster(:,indCellType) == 2,indPkTr),bigMaster(bigMaster(:,indCellType) == 2,indISI),'g.')
xlabel('Peak Trough (ms)')
ylabel('ISI Coefficient of Variation')
% title(fileName,'fontweight','bold', 'Interpreter', 'none');

%plot out proportion of each cell type

subplot(3,4,9)
cellDist = [length(findMSNs),length(findPVs),length(findCHATs)];
pie(cellDist)
labels = {'MSNs','PVs','ChATs'};
detZero = find(cellDist == 0);
labels(detZero) = [];
legend(labels,'Location','southoutside','Orientation','horizontal')



%% Column 2

subplot(3,4,2)
hold on
holder = intFastPos(findMSNs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

%plot distribution of positive, negative, both, and untuned units
subplot(3,4,6)
holder = bigMaster(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
det = holder(:,1) + holder(:,2);
det = hist(det,[-2:1:1]);
pie(det)
labels = {'Neg','Mix','None','Pos'};
detZero = find(det == 0);
labels(detZero) = [];
legend(labels,'Location','southoutside','Orientation','horizontal')
title(strcat('MSNs n=',num2str(length(findMSNs))))

subplot(3,4,10)
hold on
holder = intSlowPos(findMSNs);
hist(holder(holder > 0),widthHistVect)
xlim([0 widthHistVect(end)])
title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))
%% Column 3 PV CELLS
if findPVs
    
    subplot(3,4,4)
    hold on
    holder = intFastPos(findPVs);
    hist(holder(holder > 0),widthHistVect)
    xlim([0 widthHistVect(end)])
    title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

    %plot distribution of positive, negative, both, and untuned units
    subplot(3,4,8)
    holder = bigMaster(findPVs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    det = holder(:,1) + holder(:,2);
    det = hist(det,[-2:1:1]);
    pie(det)
    labels = {'Neg','Mix','None','Pos'};
    detZero = find(det == 0);
    labels(detZero) = [];
    legend(labels,'Location','southoutside','Orientation','horizontal')
    title(strcat('PVs n=',num2str(length(findPVs))))
    
    subplot(3,4,12)
    hold on
    holder = intSlowPos(findPVs);
    hist(holder(holder > 0),widthHistVect)
    xlim([0 widthHistVect(end)])
    title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))
end

%% Column 4 CHATs
if findCHATs
    
    subplot(3,4,3)
    hold on
    holder = intFastPos(findCHATs);
    hist(holder(holder > 0),widthHistVect)
    xlim([0 widthHistVect(end)])
    title(strcat('"int" fast pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))

    %plot distribution of positive, negative, both, and untuned units
    subplot(3,4,7)
    holder = bigMaster(findCHATs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    det = holder(:,1) + holder(:,2);
    det = hist(det,[-2:1:1]);
    pie(det)
    labels = {'Neg','Mix','None','Pos'};
    detZero = find(det == 0);
    labels(detZero) = [];
    legend(labels,'Location','southoutside','Orientation','horizontal')
    title(strcat('CHATs n=',num2str(length(findCHATs))))
    
    
    subplot(3,4,11)
    hold on
    holder = intSlowPos(findCHATs);
    hist(holder(holder > 0),widthHistVect)
    xlim([0 widthHistVect(end)])
    title(strcat('"int" slow pos, mean =',num2str(mean(holder(holder>0))),',n =',num2str(length(holder(holder>0)))))
end

spikeGraphName = 'WrapperFigure1';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



