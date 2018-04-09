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
    
    
    bigMasterInd = bigMasterInd + numUnits;
end



%now plot things!
[indCellType] = functionCellStringFind(masterHeader,'CellType');


%determine if there are units in each category
findPVs = find(bigMaster(:,indCellType) == 1);

findMSNs = find(bigMaster(:,indCellType) == 0);

findCHATs = find(bigMaster(:,indCellType) == 2);

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



