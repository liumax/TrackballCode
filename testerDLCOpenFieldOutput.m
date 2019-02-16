%%This is test code to look at DLC output from open field laser testing


cd Z:\Max\openFieldVideosSecondSet

%pull hardware files
namesHW = what;
namesHW = namesHW.mat;
hwInd = strfind(namesHW,'Hardware');
hwInd = find(not(cellfun('isempty', hwInd)));
namesHW = namesHW(hwInd);
for i = 1:length(namesHW)
    testName = namesHW{i};
    findTrial = strfind(testName,'l');
    findA = strfind(testName,'A');
    findPer = strfind(testName,'.');
    HWfinder(i,1) = str2num(testName(findTrial+1:findA-1));
    HWfinder(i,2) = str2num(testName(findPer-1));
end

%pull csv file names
namesCSV = dir;
namesCSV = {namesCSV.name}';
csvInd = strfind(namesCSV,'csv');
csvInd = find(not(cellfun('isempty', csvInd)));
namesCSV = namesCSV(csvInd);
for i = 1:length(namesCSV)
    testName = namesCSV{i};
    findTrial = strfind(testName,'Trial');
    findA = strfind(testName,'Arena');
    findPer = strfind(testName,'DeepCut');
    CSVfinder(i,1) = str2num(testName(findTrial+5:findA-1));
    CSVfinder(i,2) = str2num(testName(findPer-1));
end

%now go through and extract angle change from each
for i = 1:length(CSVfinder)
    disp(i)
    %first, lets find appropriate hardware file.
    tarFind = find(HWfinder(:,1) == CSVfinder(i,1) & HWfinder(:,2) == CSVfinder(i,2));
    dataOut = csvread(namesCSV{i},3);
    load(namesHW{tarFind})
    %calculate between the ear point.
    earMeanX = mean([dataOut(:,5),dataOut(:,8)]');
    earMeanY = mean([dataOut(:,6),dataOut(:,9)]');
    %pull angle
    for j = 1:length(earMeanY)
        v1=[earMeanX(j),earMeanY(j)]-[dataOut(j,2),dataOut(j,3)];
        v2=[earMeanX(j),earMeanY(j)]-[dataOut(j,11),dataOut(j,12)];
    %     a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
    %     angleVal(i) = a1;
        pointDir(j) = (dataOut(j,2) - earMeanX(j))*(dataOut(j,12)-earMeanY(j)) - (dataOut(j,3) - earMeanY(j))*(dataOut(j,11)-earMeanX(j));
        %positive values mean nose is left of body, negative values mean nose
        %is right of body. 
        angleVal(j)=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
    end
    %adjust to degrees
    angleVal = angleVal/(2*pi)*360;
    %correct for direction. Positive means nose is to left, negative means
    %nose is to right. 
    trueAngle = angleVal - 180;
    posNeg = zeros(length(pointDir),1);
    posNeg(pointDir > 0) = 1;
    posNeg(pointDir < 0) = -1;
    trueAngle = trueAngle.*posNeg';

    %lets pull time points of laser onset
    laserOn = timeStore(1:2:end);
    laserOn = round(laserOn*30);

    %check periods of laser on vs off
    laserDur = 10*30;
    angleLaser = [];
    angleControl = [];
    %pull cumulative angle
    for j = 1:length(laserOn)
        angleLaser(:,j) = cumsum(trueAngle(laserOn(j)-laserDur:laserOn(j) + laserDur*2));
        if laserOn(j) > laserDur*3
            angleControl(:,j) = cumsum(trueAngle(laserOn(j)-laserDur*3:laserOn(j)));
        end
    end
    
    %store these results.
    meanCumAngleLaser(:,i) = mean(angleLaser');
    STECumAngleLaser(:,i) = std(angleLaser')/sqrt(45);
    meanCumAngleControl(:,i) = mean(angleControl');
    STECumAngleLaser(:,i) = std(angleControl')/sqrt(45);
    
    angleLaser= [];
    angleControl = [];
    %pull diff angle
    for j = 1:length(laserOn)
        angleLaser(:,j) = ((trueAngle(laserOn(j)-laserDur:laserOn(j) + laserDur*2)));
        if laserOn(j) > laserDur*3
            angleControl(:,j) = ((trueAngle(laserOn(j)-laserDur*3:laserOn(j))));
        end
    end

    %store these results.
    meanDiffAngleLaser(:,i) = mean(angleLaser');
    STEDiffAngleLaser(:,i) = std(angleLaser')/sqrt(45);
    meanDiffAngleControl(:,i) = mean(angleControl');
    STEDiffAngleControl(:,i) = std(angleControl')/sqrt(45);
    
    angleLaser = [];
    angleControl = [];
end

decoder = [-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1];
% figure;
% hold on
% for i = 1:23
%     plot(meanCumAngleLaser(:,i)*decoder(i))
% end

corrCumAngleLaser = [];
for i = 1:23
    corrCumAngleLaser(:,i) = meanCumAngleLaser(:,i)*decoder(i);
end

% figure
% plot(mean(corrCumAngleLaser'))

corrDiffAngleLaser = [];
for i = 1:23
    corrDiffAngleLaser(:,i) = meanDiffAngleLaser(:,i)*decoder(i);
end

% figure
% plot(mean(corrDiffAngleLaser'))

corrCumAngleControl = [];
for i = 1:23
    corrCumAngleControl(:,i) = meanCumAngleControl(:,i)*decoder(i);
end

% figure
% plot(mean(corrCumAngleControl'))

corrDiffAngleControl = [];
for i = 1:23
    corrDiffAngleControl(:,i) = meanDiffAngleControl(:,i)*decoder(i);
end

figure
plot(mean(corrDiffAngleControl'))

figure
hold on
plot(smooth(mean(corrCumAngleLaser'),15),'g')
plot(smooth(mean(corrCumAngleControl'),15),'k')


figure
hold on
plot(smooth(mean(corrDiffAngleLaser'),15),'g')
plot(smooth(mean(corrDiffAngleControl'),15),'k')


%now we can plot while normalizing baseline
laserBase = mean(mean(corrDiffAngleLaser(100:200,:)'));
controlBase = mean(mean(corrDiffAngleControl(100:200,:)'));


hFig = figure
hold on
plot([0:1/30:30],smooth(mean(corrDiffAngleLaser'),15)-laserBase,'g','LineWidth',2)
% plot([0:1/30:30],smooth(mean(corrDiffAngleLaser')+std(corrDiffAngleLaser')/sqrt(23),15)-laserBase,'g')
% plot([0:1/30:30],smooth(mean(corrDiffAngleLaser')-std(corrDiffAngleLaser')/sqrt(23),15)-laserBase,'g')
plot([0:1/30:30],smooth(mean(corrDiffAngleControl'),15)-controlBase,'k','LineWidth',2)
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')+std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')-std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
plot([10 10],[-3 4],'k')
plot([20 20],[-3 4],'k')
spikeGraphName = 'LaserVSITI_MouseHeadAngle'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure
hold on
plot([0:1/30:30],smooth(mean(corrDiffAngleLaser'),15),'g','LineWidth',2)
plot([0:1/30:30],smooth(mean(corrDiffAngleLaser')+std(corrDiffAngleLaser')/sqrt(23),15),'g')
plot([0:1/30:30],smooth(mean(corrDiffAngleLaser')-std(corrDiffAngleLaser')/sqrt(23),15),'g')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl'),15)-controlBase,'k','LineWidth',2)
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')+std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')-std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
plot([10 10],[-3 4],'k')
plot([20 20],[-3 4],'k')
spikeGraphName = 'LaserWithSTE_MouseHeadAngle'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%% Calculate percentiles from ITI, plot these. 
%calculating percentiles by averaging across each time point, then
%performing percentile calculations. This is appropriate comparison to
%means that I am displaying. 
perLow = prctile(reshape(mean(corrDiffAngleControl(300:901,:)'),1,[]),5);
perHi = prctile(reshape(mean(corrDiffAngleControl(300:901,:)'),1,[]),95);

hFig = figure
hold on
plot([0:1/30:30],smooth(mean(corrDiffAngleLaser'),15),'g','LineWidth',2)
plot([0:1/30:30],smooth(mean(corrDiffAngleLaser')+std(corrDiffAngleLaser')/sqrt(23),15),'g')
plot([0:1/30:30],smooth(mean(corrDiffAngleLaser')-std(corrDiffAngleLaser')/sqrt(23),15),'g')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl'),15)-controlBase,'k','LineWidth',2)
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')+std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')-std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
plot([10 10],[-3 4],'k')
plot([20 20],[-3 4],'k')
plot([0 30],[perLow perLow],'k')
plot([0 30],[perHi perHi],'k')

spikeGraphName = 'LaserWithSTE_PRCTILE_MouseHeadAngle'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Try and pull longer window

%now go through and extract angle change from each
for i = 1:length(CSVfinder)
    disp(i)
    %first, lets find appropriate hardware file.
    tarFind = find(HWfinder(:,1) == CSVfinder(i,1) & HWfinder(:,2) == CSVfinder(i,2));
    dataOut = csvread(namesCSV{i},3);
    load(namesHW{tarFind})
    %calculate between the ear point.
    earMeanX = mean([dataOut(:,5),dataOut(:,8)]');
    earMeanY = mean([dataOut(:,6),dataOut(:,9)]');
    %pull angle
    for j = 1:length(earMeanY)
        v1=[earMeanX(j),earMeanY(j)]-[dataOut(j,2),dataOut(j,3)];
        v2=[earMeanX(j),earMeanY(j)]-[dataOut(j,11),dataOut(j,12)];
    %     a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
    %     angleVal(i) = a1;
        pointDir(j) = (dataOut(j,2) - earMeanX(j))*(dataOut(j,12)-earMeanY(j)) - (dataOut(j,3) - earMeanY(j))*(dataOut(j,11)-earMeanX(j));
        %positive values mean nose is left of body, negative values mean nose
        %is right of body. 
        angleVal(j)=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
    end
    %adjust to degrees
    angleVal = angleVal/(2*pi)*360;
    %correct for direction. Positive means nose is to left, negative means
    %nose is to right. 
    trueAngle = angleVal - 180;
    posNeg = zeros(length(pointDir),1);
    posNeg(pointDir > 0) = 1;
    posNeg(pointDir < 0) = -1;
    trueAngle = trueAngle.*posNeg';

    %lets pull time points of laser onset
    laserOn = timeStore(1:2:end);
    laserOn = round(laserOn*30);

    %check periods of laser on vs off
    laserDur = 10*30;
    angleLaser = [];
    %pull cumulative angle
    for j = 2:length(laserOn)
        angleLaser(:,j) = cumsum(trueAngle(laserOn(j)-2*laserDur:laserOn(j) + laserDur*3));
    end
    
    %store these results.
    longmeanCumAngleLaser(:,i) = mean(angleLaser');
    longSTECumAngleLaser(:,i) = std(angleLaser')/sqrt(45);
    
    angleLaser= [];
    %pull diff angle
    for j = 2:length(laserOn)
        angleLaser(:,j) = ((trueAngle(laserOn(j)-2*laserDur:laserOn(j) + laserDur*3)));
    end

    %store these results.
    longmeanDiffAngleLaser(:,i) = mean(angleLaser');
    longSTEDiffAngleLaser(:,i) = std(angleLaser')/sqrt(45);
    
    angleLaser = [];
end

longcorrCumAngleLaser = [];
for i = 1:23
    longcorrCumAngleLaser(:,i) = longmeanCumAngleLaser(:,i)*decoder(i);
end

figure
hold on
plot(smooth(mean(longcorrCumAngleLaser'),15),'g','LineWidth',2)
plot(smooth(mean(longcorrCumAngleLaser') + std(longcorrCumAngleLaser')/sqrt(23),15),'g')
plot(smooth(mean(longcorrCumAngleLaser') - std(longcorrCumAngleLaser')/sqrt(23),15),'g')

longcorrDiffAngleLaser = [];
for i = 1:23
    longcorrDiffAngleLaser(:,i) = longmeanDiffAngleLaser(:,i)*decoder(i);
end



longperLow = prctile(reshape(mean(longcorrDiffAngleLaser(1:600,:)'),1,[]),2.5);
longperHi = prctile(reshape(mean(longcorrDiffAngleLaser(1:600,:)'),1,[]),97.5);




hFig = figure
hold on
plot([0:1/30:50],smooth(mean(longcorrDiffAngleLaser'),15),'g','LineWidth',2)
plot([0:1/30:50],smooth(mean(longcorrDiffAngleLaser')+std(longcorrDiffAngleLaser')/sqrt(23),15),'g')
plot([0:1/30:50],smooth(mean(longcorrDiffAngleLaser')-std(longcorrDiffAngleLaser')/sqrt(23),15),'g')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl'),15)-controlBase,'k','LineWidth',2)
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')+std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')-std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
plot([20 20],[-3 4],'k')
plot([30 30],[-3 4],'k')
plot([0 50],[longperLow longperLow],'k')
plot([0 50],[longperHi longperHi],'k')

spikeGraphName = 'LongLaserWithSTE_PRCTILE_MouseHeadAngle'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

for i = 1:23
    newVals(:,i) = interp1([0:1/30:50],smooth(longcorrDiffAngleLaser(:,i),3),[0:0.1:50]);
end


hFig = figure
hold on
plot([0:0.1:50],smooth(mean(newVals'),11),'g','LineWidth',2)
plot([0:0.1:50],smooth(mean(newVals')+std(newVals')/sqrt(23),11),'g')
plot([0:0.1:50],smooth(mean(newVals')-std(newVals')/sqrt(23),11),'g')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl'),15)-controlBase,'k','LineWidth',2)
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')+std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
% plot([0:1/30:30],smooth(mean(corrDiffAngleControl')-std(corrDiffAngleControl')/sqrt(23),15)-controlBase,'k')
plot([20 20],[-3 4],'k')
plot([30 30],[-3 4],'k')
plot([0 50],[longperLow longperLow],'k')
plot([0 50],[longperHi longperHi],'k')
ylabel('Head Angle (degrees) with 95% Interval')
xlabel('Time (sec)')

spikeGraphName = 'LongLaserWithSTE_PRCTILE_MouseHeadAngle'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

% 
% hFig = figure
% hold on
% plot([0:0.1:50],smooth(cumsum(mean(newVals')),11),'g','LineWidth',2)


%% Ok lets try plotting out the allocentric body angle, not the egocentric body angle. 

%now go through and extract angle change from each
for i = 1:length(CSVfinder)
    disp(i)
    %first, lets find appropriate hardware file.
    tarFind = find(HWfinder(:,1) == CSVfinder(i,1) & HWfinder(:,2) == CSVfinder(i,2));
    dataOut = csvread(namesCSV{i},3);
    load(namesHW{tarFind})
    %calculate between the ear point.
    earMeanX = mean([dataOut(:,5),dataOut(:,8)]');
    earMeanY = mean([dataOut(:,6),dataOut(:,9)]');
    %generate vector array
    v2=[earMeanX',earMeanY']-[dataOut(:,11),dataOut(:,12)]; %This is tail
%     v2=[earMeanX',earMeanY']-[dataOut(:,2),dataOut(:,3)]; %This is tail
    preAngle = atan(v2(:,2)./v2(:,1))/(2*pi)*360;
    %eliminate non-linearities
    preDiff = diff(preAngle);
    findBig = find(abs(preDiff) > 100);
    for j = 1:length(findBig)
        %figure out the difference
        if preDiff(findBig(j)) > 0
            preAngle(findBig(j)+1:end) = preAngle(findBig(j)+1:end)- 180;
        elseif preDiff(findBig(j)) < 0
            preAngle(findBig(j)+1:end) = preAngle(findBig(j)+1:end)+ 180;
        end
    end

    %check periods of laser on vs off
    laserDur = 10*30;
    angleLaser = [];
    %pull cumulative angle
    for j = 2:length(laserOn)
%         angleLaser(:,j) = cumsum(preAngle(laserOn(j)-2*laserDur:laserOn(j) + laserDur*3))-preAngle(laserOn(j)-2*laserDur);
        angleLaser(:,j) = cumsum(preAngle(laserOn(j)-2*laserDur:laserOn(j) + laserDur*3))-preAngle(laserOn(j));
    end
    
    %store these results.
    longmeanCumAngleLaser(:,i) = mean(angleLaser');
    longSTECumAngleLaser(:,i) = std(angleLaser')/sqrt(45);
    
    angleLaser= [];
    %pull diff angle
    for j = 2:length(laserOn)
%         angleLaser(:,j) = preAngle(laserOn(j)-2*laserDur:laserOn(j) + laserDur*3)-preAngle(laserOn(j)-2*laserDur);
        angleLaser(:,j) = preAngle(laserOn(j)-2*laserDur:laserOn(j) + laserDur*3)-preAngle(laserOn(j));
    end

    %store these results.
    longmeanDiffAngleLaser(:,i) = mean(angleLaser');
    longSTEDiffAngleLaser(:,i) = std(angleLaser')/sqrt(45);
    
    angleLaser = [];
end

data15mw = [];
for i = 1:24
    data15mw(:,i) = longmeanDiffAngleLaser(:,i)*decoder(i);
end


data3mw = [];
for i = 25:48
    data3mw(:,i-24) = longmeanDiffAngleLaser(:,i)*decoder(i);
end
data3mw(:,[13,15]) = []; %dead mouse!


hFig = figure;
% longperLow = prctile(reshape(mean(data15mw(1:600,:)'),1,[]),2.5);
% longperHi = prctile(reshape(mean(data15mw(1:600,:)'),1,[]),97.5);
hold on
plot([0:1/30:50],mean(data15mw'),'g','LineWidth',2)
plot([0:1/30:50],mean(data3mw'),'b','LineWidth',2)

plot([0:1/30:50],mean(data15mw')+std(data15mw')/sqrt(24),'g','LineWidth',1)
plot([0:1/30:50],mean(data15mw')-std(data15mw')/sqrt(24),'g','LineWidth',1)
plot([0:1/30:50],mean(data3mw')+std(data3mw')/sqrt(22),'b','LineWidth',1)
plot([0:1/30:50],mean(data3mw')-std(data3mw')/sqrt(22),'b','LineWidth',1)
% plot([0 50],[longperLow longperLow],'k')
% plot([0 50],[longperHi longperHi],'k')
plot([20 20],[min(mean(data15mw')) max(mean(data15mw'))],'k')
plot([30 30],[min(mean(data15mw')) max(mean(data15mw'))],'k')
ylabel('Cumulative Mouse Axis Angle (degrees)')
xlabel('Time (sec)')
title('Cumulative Mouse Axis Angle, 15 mW (g) vs 3 mW (b)')

spikeGraphName = 'LongLaserAlloAngleWPercentile'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
smoothTar = 45;
longperLow = prctile(reshape(diff(mean(data15mw(1:600,:)')),1,[]),2.5);
longperHi = prctile(reshape(diff(mean(data15mw(1:600,:)')),1,[]),97.5);
hold on
plot([0:1/30:50-1/30],smooth(diff(mean(data15mw')),smoothTar)*30,'g','LineWidth',2)
plot([0:1/30:50-1/30],smooth(diff(mean(data15mw')),smoothTar)*30+smooth(mean(diff(data15mw)'),smoothTar)*30/sqrt(24),'g','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data15mw')),smoothTar)*30-smooth(mean(diff(data15mw)'),smoothTar)*30/sqrt(24),'g','LineWidth',1)

plot([0:1/30:50-1/30],smooth(diff(mean(data3mw')),smoothTar)*30,'b','LineWidth',2)
plot([0:1/30:50-1/30],smooth(diff(mean(data3mw')),smoothTar)*30+smooth(mean(diff(data3mw)'),smoothTar)*30/sqrt(22),'b','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data3mw')),smoothTar)*30-smooth(mean(diff(data3mw)'),smoothTar)*30/sqrt(22),'b','LineWidth',1)
plot([0 50],[longperLow*30 longperLow*30],'k')
plot([0 50],[longperHi*30 longperHi*30],'k')
plot([20 20],[min(diff(mean(data15mw')))*30 max(diff(mean(data15mw')))*30],'k')
plot([30 30],[min(diff(mean(data15mw')))*30 max(diff(mean(data15mw')))*30],'k')
ylabel('Mouse Axis Angle Change(degrees/sec) with 95% confidence')
xlabel('Time (sec)')
title('Mouse Bearing Change (degrees/sec), 15 mW (g) vs 3 mW (b)')

spikeGraphName = 'LongLaserAlloAngleChangeWPercentile'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
























figure;
hold on
for i = 1:length(decoder)
    plot(meanCumAngleLaser(:,i)*decoder(i))
end


load('HardwareTrial13Arena1.mat')
load('Trial13Arena1DLCOut.mat')
% dataOut = csvread(fname,3);
%column 1 is time, 2-3 is nose, 5-6 is left ear, 8-9 is right ear, 11-12 is
%tail base.
figure
hold on
plot(diff(dataOut(:,2)),'b')
plot(diff(dataOut(:,3)),'r')
plot((dataOut(:,4)-1)*50,'g')

figure
hold on
plot(diff(dataOut(:,5)),'b')
plot(diff(dataOut(:,6)),'r')
plot((dataOut(:,7)-1)*50,'g')


figure
hold on
plot(diff(dataOut(:,8)),'b')
plot(diff(dataOut(:,9)),'r')
plot((dataOut(:,10)-1)*50,'g')

figure
hold on
plot(diff(dataOut(:,11)),'b')
plot(diff(dataOut(:,12)),'r')
plot((dataOut(:,13)-1)*50,'g')

%seems like tail has the most shit, followed by nose. ears are doing well. 

%now lets find the angle?first, we need to generate the mean point between
%the ears.

earMeanX = mean([dataOut(:,5),dataOut(:,8)]');
earMeanY = mean([dataOut(:,6),dataOut(:,9)]');

% figure
% hold on
% tester = read(v,i);
% imagesc(tester(:,:,1))
% plot(dataOut(i,2),dataOut(i,3),'r.')
% plot(earMeanX(i),earMeanY(i),'g.')
% plot(dataOut(i,11),dataOut(i,12),'w.')


% v1=[earMeanX(i),earMeanY(i)]-[dataOut(i,2),dataOut(i,3)];
% v2=[earMeanX(i),earMeanY(i)]-[dataOut(i,11),dataOut(i,12)];


for i = 1:length(earMeanY)
    v1=[earMeanX(i),earMeanY(i)]-[dataOut(i,2),dataOut(i,3)];
    v2=[earMeanX(i),earMeanY(i)]-[dataOut(i,11),dataOut(i,12)];
%     a1 = mod(atan2( det([v1;v2;]) , dot(v1,v2) ), 2*pi );
%     angleVal(i) = a1;
    pointDir(i) = (dataOut(i,2) - earMeanX(i))*(dataOut(i,12)-earMeanY(i)) - (dataOut(i,3) - earMeanY(i))*(dataOut(i,11)-earMeanX(i));
    %positive values mean nose is left of body, negative values mean nose
    %is right of body. 
    angleVal(i)=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
end
angleVal = angleVal/(2*pi)*360;

trueAngle = angleVal - 180;
posNeg = zeros(length(pointDir),1);
posNeg(pointDir > 0) = 1;
posNeg(pointDir < 0) = -1;
trueAngle = trueAngle.*posNeg';

%lets pull time points of laser onset
laserOn = timeStore(1:2:end);
laserOn = round(laserOn*30);

%check periods of laser on vs off
laserDur = 10*30;

for i = 1:length(laserOn)
    angleLaser(:,i) = cumsum(trueAngle(laserOn(i)-laserDur:laserOn(i) + laserDur*2));
    if laserOn(i) > laserDur*3
        angleControl(:,i) = cumsum(trueAngle(laserOn(i)-laserDur*3:laserOn(i)));
    end
end


for i = 1:length(laserOn)
    angleLaser(:,i) = cumsum(diff(trueAngle(laserOn(i)-laserDur:laserOn(i) + laserDur*2)));
    if laserOn(i) > laserDur*3
        angleControl(:,i) = cumsum(diff(trueAngle(laserOn(i)-laserDur*3:laserOn(i))));
    end
end

figure
hold on
plot(mean(angleLaser'),'g')
plot(mean(angleLaser') + std(angleLaser')/sqrt(45),'g')
plot(mean(angleLaser') - std(angleLaser')/sqrt(45),'g')
plot(mean(angleControl'),'k')
plot(mean(angleControl') + std(angleControl')/sqrt(45),'k')
plot(mean(angleControl') - std(angleControl')/sqrt(45),'k')

%the readout angle doesnt seem to produce much. what about movement of the
%ear mean point?




