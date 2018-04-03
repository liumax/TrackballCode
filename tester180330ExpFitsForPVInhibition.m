%% This is code to do exponential fits for stuff from the paper. 

%% For looking at PV 2A Halo

cd Y:\Max\170905PVAnalysis
homeFolder = pwd;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
targetFolders = {'PV2AHalo/GroupData'};%,'PVARBRHalo/GroupData'};
targetBin = 0.005;
for m = 1
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
            histVect = [s.Parameters.RasterWindow(1):targetBin:s.Parameters.RasterWindow(2)];
            zeroPoint = find(histVect < 0,1,'last');
            histStore(bigInd,:) = hist(s.(s.DesignationName{j}).RasterLaser(:,1),histVect)/targetBin/length(s.Timing.LaserTimes);
%             histStore(bigInd,:) = s.(s.DesignationName{j}).HistogramLaser;
%             zWholeStore(bigInd,:) = zscore(s.(s.DesignationName{j}).HistogramLaser);
            if std(histStore(bigInd,1:zeroPoint)) > 0
                zSelfStore(bigInd,:) = (histStore(bigInd,:)-mean(histStore(bigInd,1:zeroPoint)))/std(histStore(bigInd,1:zeroPoint));
            else
                zSelfStore(bigInd,:) = zeros(length(s.(s.DesignationName{j}).HistogramLaser),1);
            end
            if mean(histStore(bigInd,1:zeroPoint)) > 0
                normSelfStore(bigInd,:) = (histStore(bigInd,:)/mean(histStore(bigInd,1:zeroPoint)));
            else
                normSelfStore(bigInd,:) = NaN(length(s.(s.DesignationName{j}).HistogramLaser),1);
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
    
    normMeanMSN = nanmean(normSelfStore(msns,:));
    normSteMSN = std(normSelfStore(msns,:))/sqrt(length(msns));
    
    pvs = find(idStore == 1);
    meanPV = mean(histStore(pvs,:));
    stePV = std(histStore(pvs,:))/sqrt(length(pvs));

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
    
end
%% Old fit system (a*e^bx + c)
% % for PV onset
% testerVect = [0:targetBin:0.05];
% tester =pvSelfZ(zeroPoint:zeroPoint+length(testerVect)-1);
% 
% ft = fittype('a*exp(b*x)+c','independent','x');
% [expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[1,-0.5,-2],...
%    'Lower',[0 -500 -6],'Upper',[20 0,0]);
% figure
% plot(expFit,testerVect,tester)
% pvOnsetFit = expFit;
% expFit
% 
% %offset. 
% 
% testerVect = [0:targetBin:0.1];
% tester =pvSelfZ(2*zeroPoint:2*zeroPoint+length(testerVect)-1);
% ft = fittype('a*exp(b*x)+c','independent','x');
% [expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-5,-0.5,-2],...
%    'Lower',[-20 -500 -2],'Upper',[0 0,2]);
% figure
% plot(expFit,testerVect,tester)
% expFit
% 
% %for MSNs onset
% 
% testerVect = [0:targetBin:0.1];
% tester =msnSelfZ(zeroPoint:zeroPoint+length(testerVect)-1);
% ft = fittype('a*exp(b*x)+c','independent','x');
% [expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-5,-0.5,0],...
%    'Lower',[-20 -500 -2],'Upper',[0 0,2]);
% figure
% plot(expFit,testerVect,tester)
% expFit
% msnOnsetFit = expFit;
% 
% 
% %offset. 
% 
% testerVect = [0:targetBin:0.4];
% tester =msnSelfZ(2*zeroPoint:2*zeroPoint+length(testerVect)-1);
% ft = fittype('a*exp(b*x)+c','independent','x');
% [expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[1,-0.5,-2],...
%    'Lower',[0 -500 -6],'Upper',[20 0,0]);
% figure
% plot(expFit,testerVect,tester)
% expFit

%% New Fit system (a*e^bx - a) or otherwise selecting mean
testerVect = [0:targetBin:0.05];
tester =pvSelfZ(zeroPoint:zeroPoint+length(testerVect)-1);

ft = fittype('a*exp(b*x)-a','independent','x');
[expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[1,-0.5],...
   'Lower',[0 -500],'Upper',[20 0]);
figure
plot(expFit,testerVect,tester)
pvOnsetFit = expFit;
expFit

%offset. 
inhibMean = mean(pvSelfZ(zeroPoint:2*zeroPoint));
testerVect = [0:targetBin:0.1];
tester =pvSelfZ(2*zeroPoint:2*zeroPoint+length(testerVect)-1);
ft = fittype('a*exp(b*x)-(a+2.0283)','independent','x');
[expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-5,-0.5],...
   'Lower',[-20 -500],'Upper',[0 0]);
figure
plot(expFit,testerVect,tester)
expFit

%for MSNs onset
inhibMean = mean(msnSelfZ(1:zeroPoint));
testerVect = [0:targetBin:0.1];
tester =msnSelfZ(zeroPoint:zeroPoint+length(testerVect)-1);
ft = fittype('a*exp(b*x)-a','independent','x');
[expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-5,-0.5],...
   'Lower',[-20 -500],'Upper',[0 0]);
figure
plot(expFit,testerVect,tester)
expFit
msnOnsetFit = expFit;


%offset. 
inhibMean = mean(msnSelfZ(zeroPoint:2*zeroPoint));
testerVect = [0:targetBin:0.4];
tester =msnSelfZ(2*zeroPoint:2*zeroPoint+length(testerVect)-1);
ft = fittype('a*exp(b*x)-a+0.3477','independent','x');
[expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[1,-0.5],...
   'Lower',[0 -500],'Upper',[20 0]);
figure
plot(expFit,testerVect,tester)
expFit

%% Do fits for normalized firing rate

%for MSNs onset
inhibMean = mean(normMeanMSN(1:zeroPoint));
testerVect = [0:targetBin:0.2];
tester =normMeanMSN(zeroPoint:zeroPoint+length(testerVect)-1);
ft = fittype('a*exp(b*x)-a+1','independent','x');
[expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-5,-0.5],...
   'Lower',[-20 -500],'Upper',[0 0]);
figure
plot(expFit,testerVect,tester)
expFit
msnOnsetFit = expFit;
% 
% 
% %offset. 
% inhibMean = mean(normMeanMSN(zeroPoint:2*zeroPoint));
% testerVect = [0:targetBin:0.4];
% tester =normMeanMSN(2*zeroPoint:2*zeroPoint+length(testerVect)-1);
% ft = fittype('a*exp(b*x)-a+0.7364','independent','x');
% [expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[1,-0.5],...
%    'Lower',[0 -500],'Upper',[20 0]);
% figure
% plot(expFit,testerVect,tester)
% expFit


%% For WT high/Low power, focusing on high power. 


cd Y:\Max\180117WTLaserControls
fileNames = what;
fileNames = fileNames.mat;
counter = 1;
targetBin = 0.005;
newWindow = [-5 10];
for i = 1:length(fileNames)
    load(fileNames{i})

    %pull master sheet
    masterDesig = s.MasterDesigs;
    masterSheet(counter:counter+size(s.MasterSheet,1)-1,:) = s.MasterSheet;
    countSheet(counter:counter+size(s.MasterSheet,1)-1) = i;
    %pull laser times
    highBlocks = s.HiBlocks;
    hiTrials = [50*(highBlocks(1) - 1)+1:50*(highBlocks(1) - 1) + 50,50*(highBlocks(2) - 1)+1:50*(highBlocks(2) - 1) + 50];
    hiLaserTimes = s.Timing.LaserTimes(hiTrials);
    %pull average rates
    for j = 1:length(s.DesignationName)
        %pull new data. 
        subhistVect = [s.Parameters.RasterWindow(1):targetBin:s.Parameters.RasterWindow(2)];
        zeroPoint = find(subhistVect < 0,1,'last');
        [newRasters] = functionBasicRaster(s.(s.DesignationName{j}).SpikeTimes,hiLaserTimes,newWindow);
        newSampleHistHiZ(counter-1+j,:) = ((hist(newRasters(:,1),subhistVect)-mean(hist(newRasters(:,1),subhistVect(1:zeroPoint))))/std(hist(newRasters(:,1),subhistVect(1:zeroPoint))))/targetBin/length(hiTrials);
        if mean(hist(newRasters(:,1),subhistVect(1:zeroPoint))) > 0
            newAverageNormHistHi(counter-1+j,:) = (hist(newRasters(:,1),subhistVect)/mean(hist(newRasters(:,1),subhistVect(1:zeroPoint))))/targetBin/length(hiTrials);
        else
            newAverageNormHistHi(counter-1+j,:) = NaN(length(subhistVect),1);
        end
        %pull old stuff too
        averageHistLow(counter-1+j,:) = s.(s.DesignationName{j}).LaserHistLow;
        averageHistHi(counter-1+j,:) = s.(s.DesignationName{j}).LaserHistHi;
        if mean(s.(s.DesignationName{j}).LaserHistHi(1:100))>0
            averageNormHistHi(counter-1+j,:) = s.(s.DesignationName{j}).LaserHistHi/mean(s.(s.DesignationName{j}).LaserHistHi(1:100));
        else
            averageNormHistHi(counter-1+j,:) = NaN(length(s.(s.DesignationName{j}).LaserHistHi),1);
        end
        averageZLow(counter-1+j,:) = (s.(s.DesignationName{j}).LaserHistLow - mean(s.(s.DesignationName{j}).HistogramLaser(1:100)))/std(s.(s.DesignationName{j}).HistogramLaser(1:100));
        averageZHi(counter-1+j,:) = (s.(s.DesignationName{j}).LaserHistHi - mean(s.(s.DesignationName{j}).HistogramLaser(1:100)))/std(s.(s.DesignationName{j}).HistogramLaser(1:100));
        isiCov(counter - 1 + j) = std(diff(s.(s.DesignationName{j}).SpikeTimes))/mean(diff(s.(s.DesignationName{j}).SpikeTimes));
    end
    counter = counter + size(s.MasterSheet,1);
end

newVect = [-5+0.025:0.05:10-0.025];

pvs = find(masterSheet(:,1) == 1 & isiCov' > 1.1);
msns = find(masterSheet(:,1) == 0 & isiCov' > 1.1);

msnLowRateMean = mean(averageHistLow(msns,:));
msnHiRateMean = mean(averageHistHi(msns,:));

msnHiRateNormMean = mean(averageNormHistHi(msns,:));
msnHiRateNormSTE = std(averageNormHistHi(msns,:))/sqrt(length(msns));

msnHiRateNewNormMean = mean(newAverageNormHistHi(msns,:));
msnHiRateNewNormSTE = std(newAverageNormHistHi(msns,:))/sqrt(length(msns));

msnLowZMean = mean(averageZLow(msns,:));
msnHiZMean = mean(averageZHi(msns,:));

pvLowRateMean = mean(averageHistLow(pvs,:));
pvHiRateMean = mean(averageHistHi(pvs,:));

pvLowZMean = mean(averageZLow(pvs,:));
pvHiZMean = mean(averageZHi(pvs,:));


%% Old exponential fits 
% %for onset
% tester =msnHiZMean(100:200);
% testerVect = [0:0.05:5];
% ft = fittype('a*exp(b*x)+c','independent','x');
% [expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-5,-0.5,0],...
%    'Lower',[0 -500 -2],'Upper',[20 0,2]);
% figure
% plot(expFit,testerVect,tester)
% expFit
% highLaserOnsetFit = expFit;
% 
% %for offset. 
% tester =msnHiZMean(200:300);
% testerVect = [0:0.05:5];
% ft = fittype('a*exp(b*x)+c','independent','x');
% [expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-1,-0.5,0],...
%    'Lower',[-20 -50 -2],'Upper',[0 0,2]);
% figure
% plot(expFit,testerVect,tester)
% expFit

steHiZ = std(averageZHi(msns,:))/sqrt(length(msns));

%% New exponential fits, constraining beginning and end
%for onset
testerVect = [0:0.005:1.5];
tester =msnHiRateNewNormMean(1000:1000+length(testerVect)-1);
testerVect = [0:0.005:1.5];
ft = fittype('a*exp(b*x)-a+1','independent','x');
[expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-5,-0.5],...
   'Lower',[0 -500],'Upper',[20 0]);
figure
plot(expFit,testerVect,tester)
expFit
highLaserOnsetFit = expFit;

%for offset. 
inhibMean = mean(msnHiRateNormMean(100:200));
tester =msnHiRateNormMean(200:300);
testerVect = [0:0.05:5];
ft = fittype('a*exp(b*x)-a-0.8116','independent','x');
[expFit,expgof] = fit(testerVect',tester',ft,'StartPoint',[-1,-0.5],...
   'Lower',[-20 -50],'Upper',[0 0]);
figure
plot(expFit,testerVect,tester)
expFit

%% Plotting




%Plot raw firing rates
figure
hold on
plot(histVect,normMeanMSN,'ko')
% plot(histVect,normMeanMSN-normSteMSN,'k--')
% plot(histVect,normMeanMSN+normSteMSN,'k--')
% plot(histVect,smooth(msnSelfZ,5),'k')
plot([0:0.005:2],msnOnsetFit([0:0.005:2]),'g','LineWidth',2)

plot(subhistVect,msnHiRateNewNormMean,'o','Color',[0.5 .5 1])
% plot(newVect,msnHiRateNormMean-msnHiRateNormSTE,'m--')
% plot(newVect,msnHiRateNormMean+msnHiRateNormSTE,'m--')
plot([0:0.05:5],highLaserOnsetFit([0:0.05:5]),'b','LineWidth',2)

xlim([-0.5 1.5])


%plot everything together
plotOffsetPV = pvOnsetFit(0);
plotOffsetMSN = msnOnsetFit(0);
plotOffsetMSNWT = highLaserOnsetFit(0);

figure
hold on
plot(histVect,pvSelfZ,'b')
plot(histVect,pvSelfZ-steSelfZPV,'b--')
plot(histVect,pvSelfZ+steSelfZPV,'b--')
% plot(histVect,smooth(pvSelfZ,5),'b')
plot([0:0.005:4],pvOnsetFit([0:0.005:4]),'r','LineWidth',2)
plot(histVect,msnSelfZ,'k')
plot(histVect,msnSelfZ-steSelfZMSN,'k--')
plot(histVect,msnSelfZ+steSelfZMSN,'k--')
% plot(histVect,smooth(msnSelfZ,5),'k')
plot([0:0.005:4],msnOnsetFit([0:0.005:4]),'g','LineWidth',2)
plot(newVect,msnHiZMean,'m')
plot(newVect,msnHiZMean-steHiZ,'m--')
plot(newVect,msnHiZMean+steHiZ,'m--')
% plot(newVect,smooth(msnHiZMean,3),'m')
plot([0:0.05:5],highLaserOnsetFit([0:0.05:5]),'c','LineWidth',2)

xlim([-0.5 1.5])


figure
hold on
% plot(histVect,pvSelfZ,'b')
% plot(histVect,smooth(pvSelfZ,5),'b')
% plot([0:0.005:4],pvOnsetFit([0:0.005:4]),'r')
plot(histVect,msnSelfZ,'k')
% plot(histVect,smooth(msnSelfZ,5),'k')
plot([0:0.005:4],msnOnsetFit([0:0.005:4]),'g')
% plot(newVect,msnHiZMean,'m')
plot(newVect,smooth(msnHiZMean,3),'m')
plot([0:0.05:5],highLaserOnsetFit([0:0.05:5]),'c')


%plot relevant onsets
figure
hold on
plot([-5:targetBin:10],smooth(mean(newSampleHistHiZ),51)) %for WT high laser
plot(histVect,pvSelfZ,'b')
plot(histVect,msnSelfZ,'k')









