%This is code to look at the completed dataset as of 170905. This dataset
%excludes ML170621A, assumes that it is a mistargeted recording.


%first, open dataset

load('170913FullDatasetWithSpikeWidthCombinedControls.mat')

%lets plot out locomotion AUC score
%first, combine all AUC storage information
bigAUCStore = [];
counter = 1;
bigAUCStore(:,1) = bigData.PV2AHalo(:,1); %store PV/MSN designation
bigAUCStore(:,2) = bigData.PV2AHalo(:,3); %store average firing rate
bigAUCStore(:,3) = bigData.PV2AHalo(:,18); %store locomotion AUC score
bigAUCStore(:,4) = bigData.PV2AHalo(:,19); %store locomotion AUC significance
bigAUCStore(:,5) = 1; %store to indicate PV2A Halo Data
counter = counter + size(bigData.PV2AHalo,1);
bigAUCStore(counter:counter+size(bigData.PVARBRHalo,1)-1,1) = bigData.PVARBRHalo(:,1); %store PV/MSN designation
bigAUCStore(counter:counter+size(bigData.PVARBRHalo,1)-1,2) = bigData.PVARBRHalo(:,3); %store average firing rate
bigAUCStore(counter:counter+size(bigData.PVARBRHalo,1)-1,3) = bigData.PVARBRHalo(:,18); %store locomotion AUC score
bigAUCStore(counter:counter+size(bigData.PVARBRHalo,1)-1,4) = bigData.PVARBRHalo(:,19); %store locomotion AUC significance
bigAUCStore(counter:counter+size(bigData.PVARBRHalo,1)-1,5) = 2; %store to indicate PV2A Halo Data
counter = counter + size(bigData.PVARBRHalo,1);

bigAUCStore(counter:counter+size(bigData.Controls,1)-1,1) = bigData.Controls(:,1); %store PV/MSN designation
bigAUCStore(counter:counter+size(bigData.Controls,1)-1,2) = bigData.Controls(:,3); %store average firing rate
bigAUCStore(counter:counter+size(bigData.Controls,1)-1,3) = bigData.Controls(:,18); %store locomotion AUC score
bigAUCStore(counter:counter+size(bigData.Controls,1)-1,4) = bigData.Controls(:,19); %store locomotion AUC significance
bigAUCStore(counter:counter+size(bigData.Controls,1)-1,5) = 3; %store to indicate PV2A Halo Data
counter = counter + size(bigData.Controls,1);

%eliminate units that arent PV or MSN
findNans = find(isnan(bigAUCStore(:,1)));
bigAUCStore(findNans,:) = [];
%eliminate any units with no running data
findNans = find(isnan(bigAUCStore(:,3)));
bigAUCStore(findNans,:) = [];
aucMSNS = find(bigAUCStore(:,1) == 0);
aucPVS = find(bigAUCStore(:,1) == 1);

msnAUCsig = intersect(aucMSNS,find(bigAUCStore(:,4)==1));
msnAUCnsig = intersect(aucMSNS,find(bigAUCStore(:,4)==0));
pvAUCsig = intersect(aucPVS,find(bigAUCStore(:,4)==1));
pvAUCnsig = intersect(aucPVS,find(bigAUCStore(:,4)==0));

msnAUCsigVals = hist(bigAUCStore(msnAUCsig,3),[0.05:0.1:0.95]);
msnAUCnsigVals = hist(bigAUCStore(msnAUCnsig,3),[0.05:0.1:0.95]);
pvAUCsigVals = hist(bigAUCStore(pvAUCsig,3),[0.05:0.1:0.95]);
pvAUCnsigVals = hist(bigAUCStore(pvAUCnsig,3),[0.05:0.1:0.95]);
%lets do some stats
n1 = length(msnAUCsig); N1 = length(aucMSNS);
n2 = length(pvAUCsig); N2 = length(aucPVS);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
hist(bigAUCStore(aucMSNS,3),[0.05:0.1:0.95])
title('MSN AUC Scores')
subplot(2,1,2)
hist(bigAUCStore(aucPVS,3),[0.05:0.1:0.95])
title('PV AUC Scores')
spikeGraphName = 'LocomotionAUCScores';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
bar([0.05:0.1:0.95],[msnAUCnsigVals;msnAUCsigVals]','Stacked')
% hist(bigAUCStore(aucMSNS,3),[0.05:0.1:0.95])
title('MSN AUC Scores')
subplot(2,1,2)
bar([0.05:0.1:0.95],[pvAUCnsigVals;pvAUCsigVals]','Stacked')
title('PV AUC Scores')
spikeGraphName = 'LocomotionAUCScoresSigVals';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,1,1)
hold on
plot(bigAUCStore(msnAUCnsig,2),bigAUCStore(msnAUCnsig,3),'k.')
plot(bigAUCStore(msnAUCsig,2),bigAUCStore(msnAUCsig,3),'r.')
title('Scatter of MSN Firing Rate vs AUC, Significant in Red')

subplot(2,1,2)
hold on
plot(bigAUCStore(pvAUCnsig,2),bigAUCStore(pvAUCnsig,3),'k.')
plot(bigAUCStore(pvAUCsig,2),bigAUCStore(pvAUCsig,3),'r.')
title('Scatter of MSN Firing Rate vs AUC, Significant in Red')
spikeGraphName = 'MsnPvAucRateScatter';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')




%plot out spike width vs other metrics
storeInd = 1;
spikeWidthRateStore = [];
subFolders = fields(bigData);
for i = 1:length(subFolders)
    holder = bigData.(subFolders{i});
    spikeWidthRateStore(storeInd:storeInd - 1 + length(holder),:) = holder;
    storeInd = storeInd + length(holder);
end

figure
plot(spikeWidthRateStore(:,2),spikeWidthRateStore(:,3),'k.')
hold on
%find significant modulation
finder = find((spikeWidthRateStore(:,8) - spikeWidthRateStore(:,6))./(spikeWidthRateStore(:,8) + spikeWidthRateStore(:,6)) < -0.5);
plot(spikeWidthRateStore(finder,2),spikeWidthRateStore(finder,3),'r.')
xlabel('Spike Width (s)')
ylabel('Firing Rate (Hz)')
title('Spike Width vs FR')

figure
plot(spikeWidthRateStore(:,2)*1000,spikeWidthRateStore(:,23)*1000,'k.')
hold on
greyFind = find(spikeWidthRateStore(:,2) > 0.0004 & spikeWidthRateStore(:,2) < 0.0005);
plot(spikeWidthRateStore(greyFind,2)*1000,spikeWidthRateStore(greyFind,23)*1000,'b.')
pvFind = find(spikeWidthRateStore(:,2) < 0.0004);
plot(spikeWidthRateStore(pvFind,2)*1000,spikeWidthRateStore(pvFind,23)*1000,'r.')
xlabel('Peak Trough (ms)')
ylabel('Spike Width (ms)')
title('Peak Trough vs Spike Width')

figure
plot(spikeWidthRateStore(:,2),spikeWidthRateStore(:,23),'k.')

%lets also include code that can basically shift the distribution of PV vs
%MSN. Try shifting to 0.4 ms

lowThresh = 0.0004;
highThresh = 0.0005;

bigData.PV2AHalo(:,1) = NaN;
bigData.PV2AHalo(bigData.PV2AHalo(:,2)<lowThresh,1) = 1;
bigData.PV2AHalo(bigData.PV2AHalo(:,2)>highThresh,1) = 0;

bigData.PVARBRHalo(:,1) = NaN;
bigData.PVARBRHalo(bigData.PVARBRHalo(:,2)<lowThresh,1) = 1;
bigData.PVARBRHalo(bigData.PVARBRHalo(:,2)>highThresh,1) = 0;

bigData.Controls(:,1) = NaN;
bigData.Controls(bigData.Controls(:,2)<lowThresh,1) = 1;
bigData.Controls(bigData.Controls(:,2)>highThresh,1) = 0;

%first test: lets see how well our controls match up!

%pull the distribution of firing rates for msns and pv cells. 

msnCtrl = find(bigData.Controls(:,1)==0);
msnHalo2A = find(bigData.PV2AHalo(:,1)==0);
msnHaloARBR = find(bigData.PVARBRHalo(:,1)==0);

pvCtrl = find(bigData.Controls(:,1)==1);
pvHalo2A = find(bigData.PV2AHalo(:,1)==1);
pvHaloARBR = find(bigData.PVARBRHalo(:,1)==1);

baseMod = [0:0.01:25];

baselineMSN = hist(bigData.Controls(msnCtrl,3),baseMod);
baselineMSN2AHalo = hist(bigData.PV2AHalo(msnHalo2A,3),baseMod);
baselineMSNARBRHalo = hist(bigData.PVARBRHalo(msnHaloARBR,3),baseMod);

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
plot(baseMod,cumsum(baselineMSN)/length(msnCtrl))
hold on
plot(baseMod,cumsum(baselineMSN2AHalo/length(msnHalo2A)),'r')
plot(baseMod,cumsum(baselineMSNARBRHalo/length(msnHaloARBR)),'m')
title('Cumulative Distribution of MSN Firing Rates')
legend(strcat('Ctrl (',num2str(length(msnCtrl)),')'),strcat('2A Halo (',num2str(length(msnHalo2A)),')'),strcat('ARBR Halo (',num2str(length(msnHaloARBR)),')'),'Location','southeast')
ylim([0 1])
xlim([0 10])

spikeGraphName = 'cumDistMSNFiring';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%based on these cumulative distributions, the control for PV2A seems to be
%a bit of an outlier.... but given variability with locomotion etc, unclear
%if this has any significance. 


baselinePVCtrl = hist(bigData.Controls(pvCtrl,3),baseMod);
% baselinePVARBRCtrl = hist(bigData.PVARBRCtrl(pvCtrlARBR,3),baseMod);
baselinePV2AHalo = hist(bigData.PV2AHalo(pvHalo2A,3),baseMod);
baselinePVARBRHalo = hist(bigData.PVARBRHalo(pvHaloARBR,3),baseMod);

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
plot(baseMod,cumsum(baselinePVCtrl)/length(pvCtrl))
hold on
plot(baseMod,cumsum(baselinePV2AHalo)/length(pvHalo2A),'r')
% plot(baseMod,cumsum(baselinePVARBRCtrl/length(pvCtrlARBR)),'g')
plot(baseMod,cumsum(baselinePVARBRHalo/length(pvHaloARBR)),'m')
title('Cumulative Distribution of PV Firing Rates')
legend(strcat('Ctrl (',num2str(length(pvCtrl)),')'),strcat('2A Halo (',num2str(length(pvHalo2A)),')'),strcat('ARBR Halo (',num2str(length(pvHaloARBR)),')'),'Location','southeast')
spikeGraphName = 'cumDistPVFiring';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%This doesnt look bad. Seems to be fairly similar among different recording
%sessions. 
vectMod = [-1:0.01:1];
%now lets look at modulation indices.
modMSNCtrl = hist((bigData.Controls(msnCtrl,8)-bigData.Controls(msnCtrl,6))./(bigData.Controls(msnCtrl,8)+bigData.Controls(msnCtrl,6)),vectMod);
modMSN2AHalo = hist((bigData.PV2AHalo(msnHalo2A,8)-bigData.PV2AHalo(msnHalo2A,6))./(bigData.PV2AHalo(msnHalo2A,8)+bigData.PV2AHalo(msnHalo2A,6)),vectMod);
% modMSNARBRCtrl = hist((bigData.PVARBRCtrl(msnCtrlARBR,8)-bigData.PVARBRCtrl(msnCtrlARBR,6))./(bigData.PVARBRCtrl(msnCtrlARBR,8)+bigData.PVARBRCtrl(msnCtrlARBR,6)),vectMod);
modMSNARBRHalo = hist((bigData.PVARBRHalo(msnHaloARBR,8)-bigData.PVARBRHalo(msnHaloARBR,6))./(bigData.PVARBRHalo(msnHaloARBR,8)+bigData.PVARBRHalo(msnHaloARBR,6)),vectMod);

modPVCtrl = hist((bigData.Controls(pvCtrl,8)-bigData.Controls(pvCtrl,6))./(bigData.Controls(pvCtrl,8)+bigData.Controls(pvCtrl,6)),vectMod);
modPV2AHalo = hist((bigData.PV2AHalo(pvHalo2A,8)-bigData.PV2AHalo(pvHalo2A,6))./(bigData.PV2AHalo(pvHalo2A,8)+bigData.PV2AHalo(pvHalo2A,6)),vectMod);
% modPVARBRCtrl = hist((bigData.PVARBRCtrl(pvCtrlARBR,8)-bigData.PVARBRCtrl(pvCtrlARBR,6))./(bigData.PVARBRCtrl(pvCtrlARBR,8)+bigData.PVARBRCtrl(pvCtrlARBR,6)),vectMod);
modPVARBRHalo = hist((bigData.PVARBRHalo(pvHaloARBR,8)-bigData.PVARBRHalo(pvHaloARBR,6))./(bigData.PVARBRHalo(pvHaloARBR,8)+bigData.PVARBRHalo(pvHaloARBR,6)),vectMod);

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
plot(vectMod,cumsum(modMSNCtrl)/length(msnCtrl))
hold on
plot(vectMod,cumsum(modMSN2AHalo/length(msnHalo2A)),'r')
% plot(vectMod,cumsum(modMSNARBRCtrl/length(msnCtrlARBR)),'g')
plot(vectMod,cumsum(modMSNARBRHalo/length(msnHaloARBR)),'m')
title('Cumulative Distribution of MSN Modulation')
legend(strcat('Ctrl (',num2str(length(msnCtrl)),')'),strcat('2A Halo (',num2str(length(msnHalo2A)),')'),strcat('ARBR Halo (',num2str(length(msnHaloARBR)),')'),'Location','southeast')
ylim([0 1])
spikeGraphName = 'cumDistMSNMOD';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
plot(vectMod,cumsum(modPVCtrl)/length(pvCtrl))
hold on
plot(vectMod,cumsum(modPV2AHalo)/length(pvHalo2A),'r')
% plot(vectMod,cumsum(modPVARBRCtrl/length(pvCtrlARBR)),'g')
plot(vectMod,cumsum(modPVARBRHalo/length(pvHaloARBR)),'m')
title('Cumulative Distribution of PV Modulation')
legend(strcat('Ctrl (',num2str(length(pvCtrl)),')'),strcat('2A Halo (',num2str(length(pvHalo2A)),')'),strcat('ARBR Halo (',num2str(length(pvHaloARBR)),')'),'Location','southeast')
ylim([0 1])
spikeGraphName = 'cumDistPVMOD';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%This doesnt look bad. Seems to be fairly similar among different recording
%sessions. 
vectMod = [-1:0.1:1];
%now lets look at modulation indices.
modMSNCtrl = hist((bigData.Controls(msnCtrl,8)-bigData.Controls(msnCtrl,6))./(bigData.Controls(msnCtrl,8)+bigData.Controls(msnCtrl,6)),vectMod);
modMSN2AHalo = hist((bigData.PV2AHalo(msnHalo2A,8)-bigData.PV2AHalo(msnHalo2A,6))./(bigData.PV2AHalo(msnHalo2A,8)+bigData.PV2AHalo(msnHalo2A,6)),vectMod);
% modMSNARBRCtrl = hist((bigData.PVARBRCtrl(msnCtrlARBR,8)-bigData.PVARBRCtrl(msnCtrlARBR,6))./(bigData.PVARBRCtrl(msnCtrlARBR,8)+bigData.PVARBRCtrl(msnCtrlARBR,6)),vectMod);
modMSNARBRHalo = hist((bigData.PVARBRHalo(msnHaloARBR,8)-bigData.PVARBRHalo(msnHaloARBR,6))./(bigData.PVARBRHalo(msnHaloARBR,8)+bigData.PVARBRHalo(msnHaloARBR,6)),vectMod);

modPVCtrl = hist((bigData.Controls(pvCtrl,8)-bigData.Controls(pvCtrl,6))./(bigData.Controls(pvCtrl,8)+bigData.Controls(pvCtrl,6)),vectMod);
modPV2AHalo = hist((bigData.PV2AHalo(pvHalo2A,8)-bigData.PV2AHalo(pvHalo2A,6))./(bigData.PV2AHalo(pvHalo2A,8)+bigData.PV2AHalo(pvHalo2A,6)),vectMod);
% modPVARBRCtrl = hist((bigData.PVARBRCtrl(pvCtrlARBR,8)-bigData.PVARBRCtrl(pvCtrlARBR,6))./(bigData.PVARBRCtrl(pvCtrlARBR,8)+bigData.PVARBRCtrl(pvCtrlARBR,6)),vectMod);
modPVARBRHalo = hist((bigData.PVARBRHalo(pvHaloARBR,8)-bigData.PVARBRHalo(pvHaloARBR,6))./(bigData.PVARBRHalo(pvHaloARBR,8)+bigData.PVARBRHalo(pvHaloARBR,6)),vectMod);


hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,1,1)
bar(vectMod,[modMSNCtrl',modMSN2AHalo'])
title('Histogram Ctrl vs 2A Halo')
xlim([-1 1])
subplot(3,1,2)
bar(vectMod,[modMSNCtrl',modMSNARBRHalo'])
title('Histogram Ctrl vs ARBR Halo')
xlim([-1 1])
subplot(3,1,3)
bar(vectMod,[modMSN2AHalo',modMSNARBRHalo'])
title('Histogram Halo 2A vs ARBR')
xlim([-1 1])
spikeGraphName = 'ComparisonHistogramsMSNsMOD';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%determine the cells that have a significant modulation pre vs laser
aVal = 0.05;
sigCellsCtrlMSN = find(bigData.Controls(:,12)<aVal & bigData.Controls(:,1) == 0);
sigCells2AMSN =  find(bigData.PV2AHalo(:,12)<aVal & bigData.PV2AHalo(:,1) == 0);
sigCellsARBRMSN =  find(bigData.PVARBRHalo(:,12)<aVal & bigData.PVARBRHalo(:,1) == 0);
insigCellsCtrlMSN = find(bigData.Controls(:,12)>=aVal & bigData.Controls(:,1) == 0);
insigCells2AMSN =  find(bigData.PV2AHalo(:,12)>=aVal & bigData.PV2AHalo(:,1) == 0);
insigCellsARBRMSN =  find(bigData.PVARBRHalo(:,12)>=aVal & bigData.PVARBRHalo(:,1) == 0);

%recalculate the histograms
vectMod = [-1:0.1:1];
inSigModCtrlMSN = hist((bigData.Controls(insigCellsCtrlMSN,8)-bigData.Controls(insigCellsCtrlMSN,6))./(bigData.Controls(insigCellsCtrlMSN,8)+bigData.Controls(insigCellsCtrlMSN,6)),vectMod);
inSigMod2AMSN = hist((bigData.PV2AHalo(insigCells2AMSN,8)-bigData.PV2AHalo(insigCells2AMSN,6))./(bigData.PV2AHalo(insigCells2AMSN,8)+bigData.PV2AHalo(insigCells2AMSN,6)),vectMod);
inSigModARBRMSN = hist((bigData.PVARBRHalo(insigCellsARBRMSN,8)-bigData.PVARBRHalo(insigCellsARBRMSN,6))./(bigData.PVARBRHalo(insigCellsARBRMSN,8)+bigData.PVARBRHalo(insigCellsARBRMSN,6)),vectMod);
SigModCtrlMSN = hist((bigData.Controls(sigCellsCtrlMSN,8)-bigData.Controls(sigCellsCtrlMSN,6))./(bigData.Controls(sigCellsCtrlMSN,8)+bigData.Controls(sigCellsCtrlMSN,6)),vectMod);
SigMod2AMSN = hist((bigData.PV2AHalo(sigCells2AMSN,8)-bigData.PV2AHalo(sigCells2AMSN,6))./(bigData.PV2AHalo(sigCells2AMSN,8)+bigData.PV2AHalo(sigCells2AMSN,6)),vectMod);
SigModARBRMSN = hist((bigData.PVARBRHalo(sigCellsARBRMSN,8)-bigData.PVARBRHalo(sigCellsARBRMSN,6))./(bigData.PVARBRHalo(sigCellsARBRMSN,8)+bigData.PVARBRHalo(sigCellsARBRMSN,6)),vectMod);

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,1,1)
bar(vectMod,[SigModCtrlMSN;inSigModCtrlMSN]','stacked')
xlim([-1 1])
title('Histogram 2A Ctrl')
subplot(3,1,2)
bar(vectMod,[SigMod2AMSN;inSigMod2AMSN]','stacked')
xlim([-1 1])
title('Histogram 2A Halo')
subplot(3,1,3)
bar(vectMod,[SigModARBRMSN;inSigModARBRMSN]','stacked')
xlim([-1 1])
title('Histogram ARBR Halo')
spikeGraphName = 'HistMSNSigChangeMod';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

vectMod = [-1:0.1:1];
histModCtrlMSN = hist((bigData.Controls(msnCtrl,8)-bigData.Controls(msnCtrl,6))./(bigData.Controls(msnCtrl,8)+bigData.Controls(msnCtrl,6)),vectMod);
histMod2AMSN = hist((bigData.PV2AHalo(msnHalo2A,8)-bigData.PV2AHalo(msnHalo2A,6))./(bigData.PV2AHalo(msnHalo2A,8)+bigData.PV2AHalo(msnHalo2A,6)),vectMod);
histModARBRMSN = hist((bigData.PVARBRHalo(msnHaloARBR,8)-bigData.PVARBRHalo(msnHaloARBR,6))./(bigData.PVARBRHalo(msnHaloARBR,8)+bigData.PVARBRHalo(msnHaloARBR,6)),vectMod);

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,1,1)
bar(vectMod,histModCtrlMSN)
xlim([-1 1])
ylim([0 100])
title('Histogram 2A Ctrl')
subplot(3,1,2)
bar(vectMod,histMod2AMSN)
xlim([-1 1])
ylim([0 100])
title('Histogram 2A Halo')
subplot(3,1,3)
bar(vectMod,histModARBRMSN)
xlim([-1 1])
ylim([0 100])
title('Histogram ARBR Halo')
spikeGraphName = 'HistMSNSimpleMod';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%determine the cells that have a significant modulation pre vs laser
aVal = 0.01;
sigCellsCtrlPV = find(bigData.Controls(:,12)<aVal & bigData.Controls(:,1) == 1);
sigCells2APV =  find(bigData.PV2AHalo(:,12)<aVal & bigData.PV2AHalo(:,1) == 1);
sigCellsARBRPV =  find(bigData.PVARBRHalo(:,12)<aVal & bigData.PVARBRHalo(:,1) == 1);
insigCellsCtrlPV = find(bigData.Controls(:,12)>=aVal & bigData.Controls(:,1) == 1);
insigCells2APV =  find(bigData.PV2AHalo(:,12)>=aVal & bigData.PV2AHalo(:,1) == 1);
insigCellsARBRPV =  find(bigData.PVARBRHalo(:,12)>=aVal & bigData.PVARBRHalo(:,1) == 1);

%recalculate the histograms
vectMod = [-1:0.1:1];
inSigModCtrlPV = hist((bigData.Controls(insigCellsCtrlPV,8)-bigData.Controls(insigCellsCtrlPV,6))./(bigData.Controls(insigCellsCtrlPV,8)+bigData.Controls(insigCellsCtrlPV,6)),vectMod);
inSigMod2APV = hist((bigData.PV2AHalo(insigCells2APV,8)-bigData.PV2AHalo(insigCells2APV,6))./(bigData.PV2AHalo(insigCells2APV,8)+bigData.PV2AHalo(insigCells2APV,6)),vectMod);
inSigModARBRPV = hist((bigData.PVARBRHalo(insigCellsARBRPV,8)-bigData.PVARBRHalo(insigCellsARBRPV,6))./(bigData.PVARBRHalo(insigCellsARBRPV,8)+bigData.PVARBRHalo(insigCellsARBRPV,6)),vectMod);
SigModCtrlPV = hist((bigData.Controls(sigCellsCtrlPV,8)-bigData.Controls(sigCellsCtrlPV,6))./(bigData.Controls(sigCellsCtrlPV,8)+bigData.Controls(sigCellsCtrlPV,6)),vectMod);
SigMod2APV = hist((bigData.PV2AHalo(sigCells2APV,8)-bigData.PV2AHalo(sigCells2APV,6))./(bigData.PV2AHalo(sigCells2APV,8)+bigData.PV2AHalo(sigCells2APV,6)),vectMod);
SigModARBRPV = hist((bigData.PVARBRHalo(sigCellsARBRPV,8)-bigData.PVARBRHalo(sigCellsARBRPV,6))./(bigData.PVARBRHalo(sigCellsARBRPV,8)+bigData.PVARBRHalo(sigCellsARBRPV,6)),vectMod);

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,1,1)
bar(vectMod,[SigModCtrlPV;inSigModCtrlPV]','stacked')
xlim([-1 1])
title('Histogram 2A Ctrl')
subplot(3,1,2)
bar(vectMod,[SigMod2APV;inSigMod2APV]','stacked')
xlim([-1 1])
title('Histogram 2A Halo')
subplot(3,1,3)
bar(vectMod,[SigModARBRPV;inSigModARBRPV]','stacked')
xlim([-1 1])
title('Histogram ARBR Halo')
spikeGraphName = 'HistPVSigChangeMod';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

vectMod = [-1:0.1:1];
histModCtrlPV = hist((bigData.Controls(pvCtrl,8)-bigData.Controls(pvCtrl,6))./(bigData.Controls(pvCtrl,8)+bigData.Controls(pvCtrl,6)),vectMod);
histMod2APV = hist((bigData.PV2AHalo(pvHalo2A,8)-bigData.PV2AHalo(pvHalo2A,6))./(bigData.PV2AHalo(pvHalo2A,8)+bigData.PV2AHalo(pvHalo2A,6)),vectMod);
histModARBRPV = hist((bigData.PVARBRHalo(pvHaloARBR,8)-bigData.PVARBRHalo(pvHaloARBR,6))./(bigData.PVARBRHalo(pvHaloARBR,8)+bigData.PVARBRHalo(pvHaloARBR,6)),vectMod);

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,1,1)
bar(vectMod,histModCtrlPV)
xlim([-1 1])
title('Histogram 2A Ctrl')
subplot(3,1,2)
bar(vectMod,histMod2APV)
xlim([-1 1])
ylim([0 6])
title('Histogram 2A Halo')
subplot(3,1,3)
bar(vectMod,histModARBRPV)
xlim([-1 1])
ylim([0 6])
title('Histogram ARBR Halo')
spikeGraphName = 'HistPVSimpleMod';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



TruemodMSNCtrl = (bigData.Controls(msnCtrl,8)-bigData.Controls(msnCtrl,6))./(bigData.Controls(msnCtrl,8)+bigData.Controls(msnCtrl,6));
TruemodMSN2AHalo = (bigData.PV2AHalo(msnHalo2A,8)-bigData.PV2AHalo(msnHalo2A,6))./(bigData.PV2AHalo(msnHalo2A,8)+bigData.PV2AHalo(msnHalo2A,6));
% TruemodMSNARBRCtrl = (bigData.PVARBRCtrl(msnCtrlARBR,8)-bigData.PVARBRCtrl(msnCtrlARBR,6))./(bigData.PVARBRCtrl(msnCtrlARBR,8)+bigData.PVARBRCtrl(msnCtrlARBR,6));
TruemodMSNARBRHalo = (bigData.PVARBRHalo(msnHaloARBR,8)-bigData.PVARBRHalo(msnHaloARBR,6))./(bigData.PVARBRHalo(msnHaloARBR,8)+bigData.PVARBRHalo(msnHaloARBR,6));

TruemodPVCtrl = (bigData.Controls(pvCtrl,8)-bigData.Controls(pvCtrl,6))./(bigData.Controls(pvCtrl,8)+bigData.Controls(pvCtrl,6));
TruemodPV2AHalo = (bigData.PV2AHalo(pvHalo2A,8)-bigData.PV2AHalo(pvHalo2A,6))./(bigData.PV2AHalo(pvHalo2A,8)+bigData.PV2AHalo(pvHalo2A,6));
% TruemodPVARBRCtrl = (bigData.PVARBRCtrl(pvCtrlARBR,8)-bigData.PVARBRCtrl(pvCtrlARBR,6))./(bigData.PVARBRCtrl(pvCtrlARBR,8)+bigData.PVARBRCtrl(pvCtrlARBR,6));
TruemodPVARBRHalo = (bigData.PVARBRHalo(pvHaloARBR,8)-bigData.PVARBRHalo(pvHaloARBR,6))./(bigData.PVARBRHalo(pvHaloARBR,8)+bigData.PVARBRHalo(pvHaloARBR,6));

% modStore = cell;
modStore(:,1) = {TruemodMSNCtrl,TruemodMSN2AHalo,TruemodMSNARBRHalo};
modStore(:,2) = {TruemodPVCtrl,TruemodPV2AHalo,TruemodPVARBRHalo};

%store gamma isis
gammaStore{1,1} = bigData.Controls(msnCtrl,22);
gammaStore{2,1} = bigData.PV2AHalo(msnHalo2A,22);
gammaStore{3,1} = bigData.PVARBRHalo(msnHaloARBR,22);

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
x = gammaStore{2,1};
y = modStore{2,1};
X = [ones(length(x),1),x];
b = X\y;
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
scatter(x,y)
hold on
plot(x,yCalc2,'--')
title('%Gamma(x) vs modulation (2AMSNS)')
spikeGraphName = 'gammaSpikevsMod2A';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
x = gammaStore{3,1};
y = modStore{3,1};
X = [ones(length(x),1),x];
b = X\y;
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
scatter(x,y)
hold on
plot(x,yCalc2,'--')
title('%Gamma(x) vs modulation (ARBRMSNS)')
spikeGraphName = 'gammaSpikevsModARBR';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%plot bars of modulation along with individual data points
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
hold on
for i = 1:3
    dataHolder = modStore{i,1};
    dataHolder(:,2) = i+ (randn(length(dataHolder),1))/10;
    bar(i,mean(dataHolder(:,1)),'FaceColor','white')
    plot(dataHolder(:,2),dataHolder(:,1),'ko')
end
meanHolder = [mean(modStore{1,1}),mean(modStore{2,1}),mean(modStore{3,1})];
steStore = [std(modStore{1,1})/sqrt(length(modStore{1,1})),std(modStore{2,1})/sqrt(length(modStore{2,1})),std(modStore{3,1})/sqrt(length(modStore{3,1}))];
errorbar(meanHolder,steStore)
title('MSN Modulation Index, Individual Points and Means')
set(gca,'XTick',[1,2,3])
set(gca,'XTickLabel',{strcat('Ctrl (',num2str(length(msnCtrl)),')'),strcat('2A Halo (',num2str(length(msnHalo2A)),')'),strcat('ARBR Halo (',num2str(length(msnHaloARBR)),')')})
ylim([-1 1])
ylabel('Modulation Index')
spikeGraphName = 'barGraphMSNModSinglePoints';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
hold on
for i = 1:3
    dataHolder = modStore{i,2};
    dataHolder(:,2) = i+ (randn(length(dataHolder),1))/10;
    bar(i,mean(dataHolder(:,1)),'FaceColor','white')
    plot(dataHolder(:,2),dataHolder(:,1),'ko')
end
meanHolder = [mean(modStore{1,2}),mean(modStore{2,2}),mean(modStore{3,2})];
steStore = [std(modStore{1,2})/sqrt(length(modStore{1,2})),std(modStore{2,2})/sqrt(length(modStore{2,2})),std(modStore{3,2})/sqrt(length(modStore{3,2}))];
errorbar(meanHolder,steStore)
title('PV Modulation Index, Individual Points and Means')
set(gca,'XTick',[1,2,3])
set(gca,'XTickLabel',{strcat('Ctrl (',num2str(length(pvCtrl)),')'),strcat('2A Halo (',num2str(length(pvHalo2A)),')'),strcat('ARBR Halo (',num2str(length(pvHaloARBR)),')')})
ylim([-1 1])
ylabel('Modulation Index')
spikeGraphName = 'barGraphPVModSinglePoints';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%DO SOME STATS. FIRST SEE IF THINGS ARE DIFFERENT FROM ZERO BY SIGNRANK

for i = 1:3
    for j = 1:2
        statResult(i,j) = signrank(modStore{i,j});
    end
end

%based on this, both controls are significantly different from zero. Both
%2A and ARBR are significantly different from zero, which is a relief. The
%2A effect is definitely more significant. With few exceptions, the ARBR
%line has more significant modulation in the controls.....which is weird
%given how the data plots out. I do wonder if its because the control looks
%like a symmetric distribution thats leftshifted, while the halo is
%irregular, and that is breaking the test. 


%now lets compare halo vs ctrl. 

for i = 1:2
    haloComp2A(i) = ranksum(modStore{1,i},modStore{2,i});
    haloCompARBR(i) = ranksum(modStore{3,i},modStore{1,i});
end

%looking at ranksum testing, it seems that both are significantly different
%from controls. 

for i = 1:2
    [haloTComp2A(i),p(i)] = ttest2(modStore{1,i},modStore{2,i});
    
    haloTCompARBR(i) = ttest2(modStore{3,i},modStore{1,i});
end

%if we do a two sample T test, we find that all differences are
%significant. 


for i = 1:2
    haloCompGENO(i) = ranksum(modStore{2,i},modStore{3,i});
end


%lets try chisquare for PV
%first, find modulation less than -0.5
chi(1,1) = length(find(modStore{1,2} < -0.5));
chi(2,1) = length(find(modStore{2,2} < -0.5));
chi(3,1) = length(find(modStore{3,2} < -0.5));

chi(1,2) = length(modStore{1,2});
chi(2,2) = length(modStore{2,2});
chi(3,2) = length(modStore{3,2});

chiPerms = nchoosek([1,2,3],2);

for i = 1:3
    N1 = chi(chiPerms(i,1),2);
    N2 = chi(chiPerms(i,2),2);
    n1 = chi(chiPerms(i,1),1);
    n2 = chi(chiPerms(i,2),1);
    
    x1 = [repmat('a',N1,1); repmat('b',N2,1)];
    x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
    
    [tbl,chi2stat,pval] = crosstab(x1,x2)
end


%This finds by ranksum that the effect size is different between MSNs but
%not PVs. YAY? I'm a big confused about the lack of difference between PVs,
%but I can see how with smaller N and much noisier effects, it wouldnt come
%out. The MSN difference is good, suggesting there is a significant effect
%of genotype. 


%okay, now lets try looking at distance vs modulation. 

%lets pull out distances and store them!
distStore = [];
distStore{1,1} = bigData.Controls(bigData.Controls(:,1) == 0,4);
distStore{1,2} = bigData.Controls(bigData.Controls(:,1) == 1,4);
distStore{2,1} = bigData.PV2AHalo(bigData.PV2AHalo(:,1) == 0,4);
distStore{2,2} = bigData.PV2AHalo(bigData.PV2AHalo(:,1) == 1,4);
distStore{3,1} = bigData.PVARBRHalo(bigData.PVARBRHalo(:,1) == 0,4);
distStore{3,2} = bigData.PVARBRHalo(bigData.PVARBRHalo(:,1) == 1,4);
% distStore = repmat(distStore,1,4);
%now calculate linear regressions

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.03 0.05], [0.03 0.01]);
for i = 1:2
    for j = 1:3
        x = distStore{j,i};
        y = modStore{j,i};
        X = [ones(length(x),1),x];
        b = X\y;
        yCalc2 = X*b;
        Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        subplot(3,2,(j-1)*2 + i)
        scatter(x,y)
        hold on
        plot(x,yCalc2,'--')
%         xlabel('Depth(um)')
%         ylabel('Modulation Index')
        title(strcat('Data with Regression',num2str(Rsq2)))
    end
end
spikeGraphName = 'linRegDepthVsMod';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%dont see any strong correlation coefficients, which is good. Suggests
%little to no relationship between depth and modulation index.

%now lets see if there is a relationship between locomotion ROC and
%modulation index.
rocStore = [];
rocStore{1,1} = bigData.Controls(bigData.Controls(:,1) == 0,18);
rocStore{1,2} = bigData.Controls(bigData.Controls(:,1) == 1,18);
rocStore{2,1} = bigData.PV2AHalo(bigData.PV2AHalo(:,1) == 0,18);
rocStore{2,2} = bigData.PV2AHalo(bigData.PV2AHalo(:,1) == 1,18);
rocStore{3,1} = bigData.PVARBRHalo(bigData.PVARBRHalo(:,1) == 0,18);
rocStore{3,2} = bigData.PVARBRHalo(bigData.PVARBRHalo(:,1) == 1,18);
% rocStore = repmat(rocStore,1,3);

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.03 0.05], [0.03 0.01]);
for i = 1:2
    for j = 1:3
        x = rocStore{j,i};
        y = modStore{j,i};
        X = [ones(length(x),1),x];
        b = X\y;
        yCalc2 = X*b;
        Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        subplot(3,2,(j-1)*2 + i)
        scatter(x,y)
        hold on
        plot(x,yCalc2,'--')
%         xlabel('Depth(um)')
%         ylabel('Modulation Index')
        title(strcat('Data with Regression',num2str(Rsq2)))
    end
end
spikeGraphName = 'linRegROCvsMod';
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%dont see any strong regression coefficients coming out. 

%% how about making z-scored firing rate plots?

%for this I need to go through the whole datasets, pull individual files,
%and from there, pull individual units. 

homeFolder = pwd;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
targetFolders = {'PV2AmCherry/GroupData','PV2AHalo/GroupData','PVARBRmCherry/GroupData','PVARBRHalo/GroupData'};
for m = 1:4
    cd (targetFolders{m})
    targets = what;
    targets = targets.mat;
    masterIndex = strfind(targets,'ML');
    masterIndex = find(not(cellfun('isempty', masterIndex)));
    targetFiles = targets(masterIndex);
    bigInd = 1;
    histStore = [];
    zWholeStore = [];
    zSelfStore = [];
    idStore = [];
    for i = 1:length(targetFiles)
        load(targetFiles{i})
        numUnits = length(s.DesignationName);
        for j =1:numUnits
            histVect = [s.Parameters.RasterWindow(1):s.Parameters.histBin:s.Parameters.RasterWindow(2)];
            zeroPoint = find(histVect < 0,1,'last');
            histStore(bigInd,:) = s.(s.DesignationName{j}).HistogramLaser;
            zWholeStore(bigInd,:) = zscore(s.(s.DesignationName{j}).HistogramLaser);
            if std(s.(s.DesignationName{j}).HistogramLaser(1:zeroPoint)) > 0
                zSelfStore(bigInd,:) = (s.(s.DesignationName{j}).HistogramLaser-mean(s.(s.DesignationName{j}).HistogramLaser(1:zeroPoint)))/std(s.(s.DesignationName{j}).HistogramLaser(1:zeroPoint));
            else
                zSelfStore(bigInd,:) = zeros(length(s.(s.DesignationName{j}).HistogramLaser),1);
            end
            widthHold = s.MasterSheet(j,2);
            if widthHold < 4*10^(-4) & std(diff(s.(s.DesignationName{j}).SpikeTimes))/mean(diff(s.(s.DesignationName{j}).SpikeTimes)) > 1.1
                idStore(bigInd) = 1;
            elseif widthHold > 5*10^(-4) & std(diff(s.(s.DesignationName{j}).SpikeTimes))/mean(diff(s.(s.DesignationName{j}).SpikeTimes)) > 1.1
                idStore(bigInd) = 0;
            else
                idStore(bigInd) = NaN;
            end
%             idStore(bigInd) = s.MasterSheet(j,1);
            bigInd = bigInd + 1;
        end
    end

    msns = find(idStore == 0);
    meanMSN = mean(histStore(msns,:));
    steMSN = std(histStore(msns,:))/sqrt(length(msns));
    
    pvs = find(idStore == 1);
    meanPV = mean(histStore(pvs,:));
    stePV = std(histStore(pvs,:))/sqrt(length(pvs));

    msnWholeZ = mean(zWholeStore(msns,:));
    pvWholeZ = mean(zWholeStore(pvs,:));

    msnSelfZ = mean(zSelfStore(msns,:));
    steSelfZMSN = std(zSelfStore(msns,:))/sqrt(length(msns));
    pvSelfZ = mean(zSelfStore(pvs,:));
    steSelfZPV = std(zSelfStore(pvs,:))/sqrt(length(pvs));
    
    zStore{m,1} = zSelfStore(msns,:);
    zStore{m,2} = zSelfStore(pvs,:);
    
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    subplot(2,2,1)
    hold on
    plot(meanMSN,'LineWidth',2)
    plot(meanMSN - steMSN)
    plot(meanMSN + steMSN)
    
    
    subplot(2,2,2)
    hold on
    plot(meanPV,'r','LineWidth',2)
    plot(meanPV - stePV)
    plot(meanPV + stePV)
    title(strcat('Averaged Firing Rates',targetFolders{m}))

    subplot(2,2,3)
    hold on
    plot(msnSelfZ,'LineWidth',2)
    plot(msnSelfZ - steSelfZMSN)
    plot(msnSelfZ + steSelfZMSN)
    ylim([-0.5 1.5])
    
    subplot(2,2,4)
    hold on
    plot(pvSelfZ,'r','LineWidth',2)
    plot(pvSelfZ - steSelfZPV)
    plot(pvSelfZ + steSelfZPV)
    title('Averaged Pre Z Rates')
    ylim([-6 1])
    
    cd (homeFolder)
    spikeGraphName = strcat('zFireRates',num2str(m));
    savefig(hFig,spikeGraphName);
    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end


%now combine controls
combCtrls{1} = [zStore{1,1};zStore{3,1}];
combCtrls{2} = [zStore{1,2};zStore{3,2}];


hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,2,1)
hold on
plot([-1.95:0.05:4],mean(combCtrls{1}),'LineWidth',2)
plot([-1.95:0.05:4],mean(combCtrls{1}) - std(combCtrls{1})/sqrt(size(combCtrls{1},1)))
plot([-1.95:0.05:4],mean(combCtrls{1}) + std(combCtrls{1})/sqrt(size(combCtrls{1},1)))
ylim([-0.5 1.5])

subplot(3,2,3)
hold on
plot([-1.95:0.05:4],mean(zStore{2,1}),'LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{2,1}) - std(zStore{2,1})/sqrt(size(zStore{2,1},1)))
plot([-1.95:0.05:4],mean(zStore{2,1}) + std(zStore{2,1})/sqrt(size(zStore{2,1},1)))
ylim([-0.5 1.5])

subplot(3,2,5)
hold on
plot([-1.95:0.05:4],mean(zStore{4,1}),'LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{4,1}) - std(zStore{4,1})/sqrt(size(zStore{4,1},1)))
plot([-1.95:0.05:4],mean(zStore{4,1}) + std(zStore{4,1})/sqrt(size(zStore{4,1},1)))
ylim([-0.5 1.5])


subplot(3,2,2)
hold on
plot([-1.95:0.05:4],mean(combCtrls{2}),'r','LineWidth',2)
plot([-1.95:0.05:4],mean(combCtrls{2}) - std(combCtrls{2})/sqrt(size(combCtrls{2},1)),'r')
plot([-1.95:0.05:4],mean(combCtrls{2}) + std(combCtrls{2})/sqrt(size(combCtrls{2},1)),'r')
ylim([-6 1])

subplot(3,2,4)
hold on
plot([-1.95:0.05:4],mean(zStore{2,2}),'r','LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{2,2}) - std(zStore{2,2})/sqrt(size(zStore{2,2},1)),'r')
plot([-1.95:0.05:4],mean(zStore{2,2}) + std(zStore{2,2})/sqrt(size(zStore{2,2},1)),'r')
ylim([-6 1])

subplot(3,2,6)
hold on
plot([-1.95:0.05:4],mean(zStore{4,2}),'r','LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{4,2}) - std(zStore{4,2})/sqrt(size(zStore{4,2},1)),'r')
plot([-1.95:0.05:4],mean(zStore{4,2}) + std(zStore{4,2})/sqrt(size(zStore{4,2},1)),'r')
ylim([-6 1])

spikeGraphName = strcat('combinedZFireRates',num2str(m));
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%now lets make some combined figures with histograms and average z-scored
%rates. 

%looks like PreZ is the best method for plotting out this information. From
%this, it becomes obvious that there is no effect in the controls, and is a
%stronger effect in the 2A vs the ARBR. 


%now lets try and pull waveforms from all datasheets. This way we can pick
%out a couple of example sets or make average waveforms for 

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,2,1)
hold on
plot([-1.95:0.05:4],mean(combCtrls{1}),'LineWidth',2)
plot([-1.95:0.05:4],mean(combCtrls{1}) - std(combCtrls{1})/sqrt(size(combCtrls{1},1)))
plot([-1.95:0.05:4],mean(combCtrls{1}) + std(combCtrls{1})/sqrt(size(combCtrls{1},1)))
ylim([-0.5 1.5])

subplot(3,2,3)
hold on
plot([-1.95:0.05:4],mean(zStore{2,1}),'LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{2,1}) - std(zStore{2,1})/sqrt(size(zStore{2,1},1)))
plot([-1.95:0.05:4],mean(zStore{2,1}) + std(zStore{2,1})/sqrt(size(zStore{2,1},1)))
ylim([-0.5 1.5])

subplot(3,2,5)
hold on
plot([-1.95:0.05:4],mean(zStore{4,1}),'LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{4,1}) - std(zStore{4,1})/sqrt(size(zStore{4,1},1)))
plot([-1.95:0.05:4],mean(zStore{4,1}) + std(zStore{4,1})/sqrt(size(zStore{4,1},1)))
ylim([-0.5 1.5])


subplot(3,2,2)
hold on
bar(vectMod,histModCtrlMSN)
xlim([-1 1])
title('Histogram 2A Ctrl')

subplot(3,2,4)
hold on
bar(vectMod,histMod2AMSN)
xlim([-1 1])
title('Histogram 2A Halo')

subplot(3,2,6)
bar(vectMod,histModARBRMSN)
xlim([-1 1])
title('Histogram ARBR Halo')

spikeGraphName = strcat('combinedZFireRateAndHistMSN',num2str(m));
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
subplot(3,2,1)
hold on
plot([-1.95:0.05:4],mean(combCtrls{2}),'r','LineWidth',2)
plot([-1.95:0.05:4],mean(combCtrls{2}) - std(combCtrls{2})/sqrt(size(combCtrls{2},1)),'r')
plot([-1.95:0.05:4],mean(combCtrls{2}) + std(combCtrls{2})/sqrt(size(combCtrls{2},1)),'r')
ylim([-6 1])

subplot(3,2,3)
hold on
plot([-1.95:0.05:4],mean(zStore{2,2}),'r','LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{2,2}) - std(zStore{2,2})/sqrt(size(zStore{2,2},1)),'r')
plot([-1.95:0.05:4],mean(zStore{2,2}) + std(zStore{2,2})/sqrt(size(zStore{2,2},1)),'r')
ylim([-6 1])

subplot(3,2,5)
hold on
plot([-1.95:0.05:4],mean(zStore{4,2}),'r','LineWidth',2)
plot([-1.95:0.05:4],mean(zStore{4,2}) - std(zStore{4,2})/sqrt(size(zStore{4,2},1)),'r')
plot([-1.95:0.05:4],mean(zStore{4,2}) + std(zStore{4,2})/sqrt(size(zStore{4,2},1)),'r')
ylim([-6 1])


subplot(3,2,2)
hold on
bar(vectMod,histModCtrlPV)
xlim([-1 1])
title('Histogram 2A Ctrl')

subplot(3,2,4)
hold on
bar(vectMod,histMod2APV)
xlim([-1 1])
title('Histogram 2A Halo')

subplot(3,2,6)
hold on
bar(vectMod,histModARBRPV)
xlim([-1 1])
title('Histogram ARBR Halo')

spikeGraphName = strcat('combinedZFireRateAndHistPV',num2str(m));
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')







% %lets take another look at velocity, to determine whether or not there are
% %subtle changes to the velocity that may explain the change in control
% %behavior.
% 
% 
% homeFolder = pwd;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
% targetFolders = {'PV2AmCherry/GroupData','PV2AHalo/GroupData','PVARBRmCherry/GroupData','PVARBRHalo/GroupData'};
% for m = 1:4
%     cd (targetFolders{m})
%     targets = what;
%     targets = targets.mat;
%     masterIndex = strfind(targets,'ML');
%     masterIndex = find(not(cellfun('isempty', masterIndex)));
%     targetFiles = targets(masterIndex);
%     bigInd = 1;
%     velRaster = [];
%     for i = 1:length(targetFiles)
%         load(targetFiles{i})
%         velTrace = s.RotaryData.Velocity;
%         alignTimes = s.LaserData.LaserStartTimes;
%         for j = 1:length(alignTimes)
%             matchTime = find(velTrace(:,1) - alignTimes(j) > 0,1,'first');
%             if matchTime
%                 velRaster(bigInd,:) = velTrace(matchTime - 200:matchTime + 400,2);
%                 bigInd = bigInd + 1;
%             end
%         end
%     end
%     
%     hFig = figure;
%     set(hFig, 'Position', [10 80 1240 850])
%     hold on
%     plot(mean(velRaster),'LineWidth',2)
%     plot(mean(velRaster) - std(velRaster)/sqrt(length(velRaster)))
%     plot(mean(velRaster) + std(velRaster)/sqrt(length(velRaster)))
%     xlim([1,601])
%     set(gca,'XTick',[1:100:601])
%     set(gca,'XTickLabel',[-2:1:4])
%     title(targetFolders{m})
%     
%     cd (homeFolder)
%     spikeGraphName = strcat('combVelTrace',num2str(m));
%     savefig(hFig,spikeGraphName);
%     %save as PDF with correct name
%     set(hFig,'Units','Inches');
%     pos = get(hFig,'Position');
%     set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(hFig,spikeGraphName,'-dpdf','-r0')
% end
% 
% %lol dont see any coherent change 
% 


