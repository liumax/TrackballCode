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


decoder = [-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,...
    -1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1];

%% Look through and try and find body angle. 
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

decoder = [-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,...
    -1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1];
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

% Calculate percentiles from ITI, plot these. 
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

data15mw = [];
for i = 1:24
    data15mw(:,i) = longmeanDiffAngleLaser(:,i)*decoder(i);
end


data3mw = [];
for i = 25:48
    data3mw(:,i-24) = longmeanDiffAngleLaser(:,i)*decoder(i);
end
data3mw(:,[13,15]) = []; %dead mouse!

data15mwPulse = [];
for i = 49:72
    data15mwPulse(:,i-48) = longmeanDiffAngleLaser(:,i)*decoder(i);
end
data15mwPulse(:,[13,15]) = []; %dead mouse!

data67mw = [];
for i = 73:96
    data67mw(:,i-72) = longmeanDiffAngleLaser(:,i)*decoder(i);
end
data67mw(:,[13,15]) = []; %dead mouse!


hFig = figure;
% longperLow = prctile(reshape(mean(data15mw(1:600,:)'),1,[]),2.5);
% longperHi = prctile(reshape(mean(data15mw(1:600,:)'),1,[]),97.5);
hold on
smoothTar = 45;
plot([0:1/30:50],smooth(mean(data15mw'),smoothTar),'r','LineWidth',2)
plot([0:1/30:50],smooth(mean(data67mw'),smoothTar),'Color',[0.5 0 0],'LineWidth',2)
plot([0:1/30:50],smooth(mean(data3mw'),smoothTar),'k','LineWidth',2)
plot([0:1/30:50],smooth(mean(data15mwPulse'),smoothTar),'c--','LineWidth',2)


% plot([0:1/30:50],mean(data15mw')+std(data15mw')/sqrt(24),'g','LineWidth',1)
% plot([0:1/30:50],mean(data15mw')-std(data15mw')/sqrt(24),'g','LineWidth',1)
% plot([0:1/30:50],mean(data3mw')+std(data3mw')/sqrt(22),'b','LineWidth',1)
% plot([0:1/30:50],mean(data3mw')-std(data3mw')/sqrt(22),'b','LineWidth',1)
% plot([0 50],[longperLow longperLow],'k')
% plot([0 50],[longperHi longperHi],'k')
plot([20 20],[min(mean(data15mw')) max(mean(data15mw'))],'k')
plot([30 30],[min(mean(data15mw')) max(mean(data15mw'))],'k')
ylabel('Mouse Body Angle (degrees)')
xlabel('Time (sec)')
title('Mouse Body Angle, 3, 6.7, 15 mW, 15 mW pulsed')





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
%     figure
%     hist(diff(preAngle),[-30:1:30])
%     xlim([-30 30])
    threshVal = 15;
    whileTrig = 0;
    advanceCount = 1;
    velX = diff(preAngle);
    while whileTrig == 0
        if advanceCount == length(velX)
            break
        end
        %check absolute size of frontVel
        velVal = velX(advanceCount);
        if abs(velVal) > threshVal
            disp(strcat('Threshold Crossing at',num2str(advanceCount)))
            disp(strcat('Size of',num2str(velVal)))
            if velVal > 0
                disp('Positive Threshold Crossing')
                %next we need to determine if there is a return
            elseif velVal < 0
                disp('Negative Threshold Crossing')
            end
            %now lets go and crawl through the next bit of the data. Go up
            %to 100 far. 
%             absData = abs(velX(advanceCount:advanceCount + 100));
            newCount = 0;
            while whileTrig == 0
                newVal = abs(velX(advanceCount+1+newCount) - velX(advanceCount+newCount));
                if newVal > threshVal;
                    disp('STILL BIG')
                    newCount = newCount + 1;
                else
                    disp('Found END')
                    tarCount = newCount;
                    break
                end
            end
            %now lets delete the data in between that is shitty.
            preAngle(advanceCount+1:advanceCount+tarCount+1) = linspace(preAngle(advanceCount),preAngle(advanceCount+tarCount+2),tarCount+1);
            advanceCount = advanceCount + tarCount + 1;
        else
            advanceCount = advanceCount + 1;
        end
    end
    
    %now lets do some shitty ass noldus rotation tracking
    noldusSmooth = 31;
    smoothAngle = smooth(preAngle,noldusSmooth);
    whileTrig = 0;
    creepInd = 1;
    angleTracker = zeros(length(smoothAngle),1);
%     angleInd = 1;
    angleDel = diff(smoothAngle);
    angleThresh = 90;
    angleCorr = 45;
    signVal = 0;
    while whileTrig == 0
        if creepInd >= length(angleDel)
            break
        end
        %determine sign of angle change. 
        signVal = sign(angleDel(creepInd));
        %deal with potential zeros.
        if signVal == 0
            nextFind = find(angleDel(creepInd:end) ~=0,1,'first');
            creepInd = creepInd + nextFind - 1;
        end
        %find next zero crossing for derivative. 
        if signVal > 0
            nextFind = find(angleDel(creepInd:end) < 0,1,'first');
        elseif signVal < 0
            nextFind = find(angleDel(creepInd:end) > 0,1,'first');
        end
        %reindex this properly.
        nextFind = creepInd + nextFind - 1;
        if isempty(nextFind)
            break
        end
        %the start point to the end point determine the target region. 
        %now we need to determine if there is at least 90 degrees in the
        %target region. 
        angleChange = preAngle(nextFind) - preAngle(creepInd);
        if abs(angleChange) < angleThresh %this is the case where there is insufficient rotation
            disp(strcat('Insufficient rotation,creepInd:',num2str(creepInd)))
            disp(num2str(angleChange))
            if ~isempty(nextFind)
                creepInd = nextFind;
            else
                creepInd = length(preAngle);
            end
        elseif abs(angleChange) >= angleThresh 
            disp(strcat('Sufficient Rotation',num2str(angleChange)))
            %now we need to go about marking rotations. 
            diffVals = abs(preAngle(creepInd:nextFind)-preAngle(creepInd));
            while whileTrig == 0
                %find the next 90 degrees
                find90 = find(diffVals > angleThresh,1,'first');
                if find90 %if get 90 degrees.
                    disp(strcat('90 degree rotation found:',num2str(creepInd),'->',num2str(creepInd + find90)))
                    %saves data point
                    angleTracker(creepInd + find90 - 1) = signVal;
                    creepInd = creepInd + find90;
                    diffVals = abs(preAngle(creepInd:nextFind)-preAngle(creepInd));
                else %if dont have 90. 
                    %now we need to find next zero crossing of derivative
                    if angleDel(nextFind) > 0
                        nextNextFind = find(angleDel(nextFind:end) < 0,1,'first');
                    elseif angleDel(nextFind) < 0
                        nextNextFind = find(angleDel(nextFind:end) > 0,1,'first');
                    end
                    %now we need to examine the angle value accumulated
                    %over that period
                    nextAngle = preAngle(nextFind + nextNextFind - 1) - preAngle(nextFind);
                    if abs(nextAngle) > angleCorr %if this is bigger. 
                        creepInd = nextFind;
                        disp('Interruption Detected, Greater than Threshold.')
                        break %exit small while loop to restart
                    else %if less  than threshold. then we want to jump to the next zero crossing. 
                        disp('Interruption Detected, Less than Threshold, skipping over')
                        disp(strcat('CreepInd:',num2str(creepInd),'nextFind:',num2str(nextFind)))
                        nextFind = nextFind + nextNextFind -1;
                        if angleDel(nextFind) > 0
                            nextNextFind = find(angleDel(nextFind:end) < 0,1,'first');
                        elseif angleDel(nextFind) < 0
                            nextNextFind = find(angleDel(nextFind:end) > 0,1,'first');
                        end
                        if nextNextFind
                            diffVals = abs(preAngle(creepInd:nextFind+nextNextFind-1)-preAngle(creepInd));
                            disp(strcat('New Diff:',num2str(diffVals(end))))
                        else %this means you can find another zero crossing anymore. 
                            creepInd = length(angleDel);
                            break
                        end
                    end
                    
                end
                
            end
        end
    end
    
    figure
    plot(preAngle)
    hold on
    plot(cumsum(angleTracker*90),'r')
    
    %lets pull time points of laser onset
    laserOn = timeStore(1:2:end);
    laserOn = round(laserOn*30);

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
    
    %pull rotation data. 
    rotStore = [];
%     rotNEG = [];
    basePer = [-laserDur/2 0];
    laserPer = [laserDur/2 laserDur];
    postPer = [laserDur+laserDur/2 2*laserDur];
    for j = 1:length(laserOn)
        %find positive and negative rotations in baseline period
        prePos = length(find(angleTracker(laserOn(j)+basePer(1):laserOn(j)+basePer(2)) == 1));
        preNeg = length(find(angleTracker(laserOn(j)+basePer(1):laserOn(j)+basePer(2)) == -1));
        laserPos = length(find(angleTracker(laserOn(j)+laserPer(1):laserOn(j)+laserPer(2)) == 1));
        laserNeg = length(find(angleTracker(laserOn(j)+laserPer(1):laserOn(j)+laserPer(2)) == -1));
        postPos = length(find(angleTracker(laserOn(j)+postPer(1):laserOn(j)+postPer(2)) == 1));
        postNeg = length(find(angleTracker(laserOn(j)+postPer(1):laserOn(j)+postPer(2)) == -1));
        rotStore(j,:) = [prePos,preNeg,laserPos,laserNeg,postPos,postNeg];
    end
    bigRotStore(i,:) = sum(rotStore);
end

% %normalize rotations. This is contra over all
% normRot = [];
% for i = 1:length(bigRotStore)
%     if decoder(i) == -1
%         normRot(i,:) = [bigRotStore(i,2)/(bigRotStore(i,2)+bigRotStore(i,1)),bigRotStore(i,4)/(bigRotStore(i,3)+bigRotStore(i,4)),bigRotStore(i,6)/(bigRotStore(i,5)+bigRotStore(i,6))];
%     else
%         normRot(i,:) = [bigRotStore(i,1)/(bigRotStore(i,2)+bigRotStore(i,1)),bigRotStore(i,3)/(bigRotStore(i,3)+bigRotStore(i,4)),bigRotStore(i,5)/(bigRotStore(i,5)+bigRotStore(i,6))];
%     end
% end

%alternative. can do subtraction over addition
normRot = [];
for i = 1:length(bigRotStore)
    if decoder(i) == -1
        normRot(i,:) = [(bigRotStore(i,2)-bigRotStore(i,1))/(bigRotStore(i,2)+bigRotStore(i,1)),(bigRotStore(i,4)-bigRotStore(i,3))/(bigRotStore(i,3)+bigRotStore(i,4)),(bigRotStore(i,6)-bigRotStore(i,5))/(bigRotStore(i,5)+bigRotStore(i,6))];
    else
        normRot(i,:) = [(bigRotStore(i,1)-bigRotStore(i,2))/(bigRotStore(i,2)+bigRotStore(i,1)),(bigRotStore(i,3)-bigRotStore(i,4))/(bigRotStore(i,3)+bigRotStore(i,4)),(bigRotStore(i,5)-bigRotStore(i,6))/(bigRotStore(i,5)+bigRotStore(i,6))];
    end
end

%eliminate dead mouse
targets = [13,15,13+24,15+24,13+48,15+48,13+72,15+72];
normRot(targets,:) = NaN;

% normRot = normRot./normRot(:,1);
normRotNorm(:,1) = normRot(:,1)./normRot(:,1);
normRotNorm(:,2) = normRot(:,2)./normRot(:,1);
normRotNorm(:,3) = normRot(:,3)./normRot(:,1);
figure
hold on
for i = 1:24
    plot(normRotNorm(i,:))
end

rots15mW = normRot(1:24,:);
rots3mW = normRot(25:48,:);
rots15mWPulse = normRot(49:72,:);
rots7mW = normRot(73:96,:);


rots15mW([13,15],:) = [];
rots3mW([13,15],:) = [];
rots15mWPulse([13,15],:) = [];
rots7mW([13,15],:) = [];

signrank(mean([rots3mW(:,1),rots7mW(:,3)]'),rots7mW(:,2))
signrank(mean([rots7mW(:,1),rots7mW(:,3)]'),rots7mW(:,2))
signrank(mean([rots15mW(:,1),rots7mW(:,3)]'),rots7mW(:,2))
signrank(mean([rots15mWPulse(:,1),rots7mW(:,3)]'),rots7mW(:,2))

hFig = figure;
hold on
plot([1,2,3],rots3mW,'k')
errorbar([1,2,3],mean(rots3mW),std(rots3mW)/sqrt(length(rots3mW)),'r','LineWidth',2)
plot([4,5,6],rots7mW,'k')
errorbar([4,5,6],mean(rots7mW),std(rots7mW)/sqrt(length(rots7mW)),'r','LineWidth',2)
plot([7,8,9],rots15mW,'k')
errorbar([7,8,9],mean(rots15mW),std(rots15mW)/sqrt(length(rots15mW)),'r','LineWidth',2)
plot([10,11,12],rots15mWPulse,'k')
errorbar([10,11,12],mean(rots15mWPulse),std(rots15mWPulse)/sqrt(length(rots15mWPulse)),'r','LineWidth',2)
ylim([-1 1])
set(gca,'tickdir','out')

spikeGraphName = 'DLCRotationsScatterPlot'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


data15mw = [];
for i = 1:24
    data15mw(:,i) = longmeanDiffAngleLaser(:,i)*decoder(i);
end


data3mw = [];
for i = 25:48
    data3mw(:,i-24) = longmeanDiffAngleLaser(:,i)*decoder(i);
end
data3mw(:,[13,15]) = []; %dead mouse!

data15mwPulse = [];
for i = 49:72
    data15mwPulse(:,i-48) = longmeanDiffAngleLaser(:,i)*decoder(i);
end
data15mwPulse(:,[13,15]) = []; %dead mouse!

data67mw = [];
for i = 73:96
    data67mw(:,i-72) = longmeanDiffAngleLaser(:,i)*decoder(i);
end
data67mw(:,[13,15]) = []; %dead mouse!

hFig = figure;
% longperLow = prctile(reshape(mean(data15mw(1:600,:)'),1,[]),2.5);
% longperHi = prctile(reshape(mean(data15mw(1:600,:)'),1,[]),97.5);
hold on
plot([0:1/30:50],mean(data15mw'),'r','LineWidth',2)
plot([0:1/30:50],mean(data67mw'),'Color',[0.5 0 0],'LineWidth',2)
plot([0:1/30:50],mean(data3mw'),'k','LineWidth',2)
plot([0:1/30:50],mean(data15mwPulse'),'c--','LineWidth',2)


% plot([0:1/30:50],mean(data15mw')+std(data15mw')/sqrt(24),'r','LineWidth',1)
% plot([0:1/30:50],mean(data15mw')-std(data15mw')/sqrt(24),'r','LineWidth',1)
% plot([0:1/30:50],mean(data3mw')+std(data3mw')/sqrt(22),'k','LineWidth',1)
% plot([0:1/30:50],mean(data3mw')-std(data3mw')/sqrt(22),'k','LineWidth',1)
% plot([0:1/30:50],mean(data67mw')+std(data67mw')/sqrt(22),'Color',[0.5 0 0],'LineWidth',1)
% plot([0:1/30:50],mean(data67mw')-std(data67mw')/sqrt(22),'Color',[0.5 0 0],'LineWidth',1)
% plot([0:1/30:50],mean(data15mwPulse')+std(data15mwPulse')/sqrt(22),'c--','LineWidth',1)
% plot([0:1/30:50],mean(data15mwPulse')-std(data15mwPulse')/sqrt(22),'c--','LineWidth',1)

plot([20 20],[min(mean(data15mw')) max(mean(data15mw'))],'k')
plot([30 30],[min(mean(data15mw')) max(mean(data15mw'))],'k')
ylabel('Mouse Body Direction (degrees)')
xlabel('Time (sec)')
title('Mouse Body Direction, 3, 6.7, 15 mW, 15 mW pulsed')

spikeGraphName = 'LongLaserAlloAngleWPercentile'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
smoothTar = 45;
longperLow = prctile(reshape([diff(mean(data15mw(1:600,:)')),diff(mean(data3mw(1:600,:)')),diff(mean(data67mw(1:600,:)')),diff(mean(data15mwPulse(1:600,:)'))],1,[]),2.5);
longperHi = prctile(reshape([diff(mean(data15mw(1:600,:)')),diff(mean(data3mw(1:600,:)')),diff(mean(data67mw(1:600,:)')),diff(mean(data15mwPulse(1:600,:)'))],1,[]),97.5);
hold on
plot([0:1/30:50-1/30],smooth(diff(mean(data15mw')),smoothTar)*30,'r','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data15mw')),smoothTar)*30+smooth(mean(diff(data15mw)'),smoothTar)*30/sqrt(24),'r','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data15mw')),smoothTar)*30-smooth(mean(diff(data15mw)'),smoothTar)*30/sqrt(24),'r','LineWidth',1)

plot([0:1/30:50-1/30],smooth(diff(mean(data15mwPulse')),smoothTar)*30,'c','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data15mwPulse')),smoothTar)*30+smooth(mean(diff(data15mwPulse)'),smoothTar)*30/sqrt(22),'c','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data15mwPulse')),smoothTar)*30-smooth(mean(diff(data15mwPulse)'),smoothTar)*30/sqrt(22),'c','LineWidth',1)

plot([0:1/30:50-1/30],smooth(diff(mean(data3mw')),smoothTar)*30,'k','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data3mw')),smoothTar)*30+smooth(mean(diff(data3mw)'),smoothTar)*30/sqrt(22),'k','LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data3mw')),smoothTar)*30-smooth(mean(diff(data3mw)'),smoothTar)*30/sqrt(22),'k','LineWidth',1)

plot([0:1/30:50-1/30],smooth(diff(mean(data67mw')),smoothTar)*30,'Color',[0.5 0 0],'LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data67mw')),smoothTar)*30+smooth(mean(diff(data67mw)'),smoothTar)*30/sqrt(22),'Color',[0.5 0 0],'LineWidth',1)
plot([0:1/30:50-1/30],smooth(diff(mean(data67mw')),smoothTar)*30-smooth(mean(diff(data67mw)'),smoothTar)*30/sqrt(22),'Color',[0.5 0 0],'LineWidth',1)
plot([0 50],[longperLow*30 longperLow*30],'k')
plot([0 50],[longperHi*30 longperHi*30],'k')
plot([20 20],[min(diff(mean(data15mw')))*30 max(diff(mean(data15mw')))*30],'k')
plot([30 30],[min(diff(mean(data15mw')))*30 max(diff(mean(data15mw')))*30],'k')
ylabel('Mouse Axis Angle Change(degrees/sec) with 95% confidence')
xlabel('Time (sec)')
title('Mouse Bearing Change (degrees/sec), 3 mW, 6. mW, 15 mW, 15 mWPulsed')

spikeGraphName = 'LongLaserAlloAngleChangeWPercentile'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%generate figure of tracking
v = VideoReader('Trial13Arena1.mp4')
load('HardwareTrial13Arena1.mat')
load('Trial13Arena1DLCOut.mat')
tester = read(v,1);
tester = mean(tester,3);

hFig = figure;
imagesc(tester)
colormap('gray')
hold on
p1 = [dataOut(1,11) dataOut(1,12)];                         % First Point
p2 = [mean([dataOut(1,5),dataOut(1,8)]) mean([dataOut(1,6),dataOut(1,9)])]; % Second Point
dp = p2-p1;
% plot(p1,'r.')
% plot(p2,'g.')
quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'Color',[1 0 0])
spikeGraphName = 'SingleFrameAxisExample'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
imagesc(tester)
colormap('gray')
hold on
frameLim = 1;
frameMult = 1; %this instructs whether to jump frames. 
for i = 1:frameLim
    p1 = [dataOut(i*frameMult,11) dataOut(i*frameMult,12)];                         % First Point
    p2 = [mean([dataOut(i*frameMult,5),dataOut(i*frameMult,8)]) mean([dataOut(i*frameMult,6),dataOut(i*frameMult,9)])]; % Second Point
    dp = p2-p1;
    % plot(p1,'r.')
    % plot(p2,'g.')
    quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'Color',[(frameLim-i)/frameLim i/frameLim 0])
end
spikeGraphName = 'SingleFrameAxisExample'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%set parameters.
frameLim = 30;
frameMult = 2; %this instructs whether to jump frames. 
%find last frame.
newFrame = read(v,frameLim*frameMult);
newFrame = mean(newFrame,3);

combIm = (tester+newFrame)/2;

hFig = figure;
imagesc(combIm)
colormap('gray')
hold on

for i = 1:frameLim
    p1 = [dataOut(i*frameMult,11) dataOut(i*frameMult,12)];                         % First Point
    p2 = [mean([dataOut(i*frameMult,5),dataOut(i*frameMult,8)]) mean([dataOut(i*frameMult,6),dataOut(i*frameMult,9)])]; % Second Point
    dp = p2-p1;
    % plot(p1,'r.')
    % plot(p2,'g.')
    quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'Color',[(frameLim-i)/frameLim i/frameLim 0])
end
spikeGraphName = 'MultiFrameAxisExample'
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



plot([],[],'r')

p1 = [2 3];                         % First Point
p2 = [9 8];                         % Second Point
dp = p2-p1;                         % Difference
figure(1)
quiver(p1(1),p1(2),dp(1),dp(2),0)




















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




