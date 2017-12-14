%This is meant to be wrapper function for analysis of pavlovian behavior
%with photometry. This should be executed on a per animal basis. 

clear

%set parameters
behavLim = 0.8;

targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);
masterIndex = strfind(targetFiles,'Analysis');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);


numFiles = length(targetFiles);
incInd = 1;
for bigInd = 1:numFiles
    bigInd
    
    %load the target file
    load(targetFiles{bigInd})
    
    %first thing we want to do is pull out all the fluorescence traces and
    %pop them into a giant array. We also want to pull out the assumed peak
    %value from each of these, so that we have a more simple metric to
    %compare. I want to pull peaks from onset and offset. 
    
    photoTimeStep = mean(diff(s.Photo.Photo.x70dFTime));
    rasterVect = [-3:photoTimeStep:5];
    zeroPoint = find(rasterVect > 0,1,'first');
    endPoint = find(rasterVect > 1,1,'first');
    
    preBins = reshape(s.PhotoRaster.ToneRaster([1:zeroPoint],:),1,[]);
    preMean = mean(preBins);
    preSTD = std(preBins);
    zRaster = (s.PhotoRaster.ToneRaster - preMean)/preSTD;
    
    %baseline period will simply be average of the five bins before zero
    zeroBase = mean(zRaster(zeroPoint - 4:zeroPoint,:));
    endBase =  mean(zRaster(endPoint - 4:endPoint,:));
    
    %find peak values for the twentyfive bins following the stimulus onset
    %(~500ms)
    zeroPeak = max(zRaster(zeroPoint:zeroPoint+ 25,:));
    endPeak = max(zRaster(endPoint:endPoint+ 25,:));
    
    %capture average velocity from different time points
    preVel = mean(s.VelRaster.ToneRaster([19:29],:));
    zeroVel = mean(s.VelRaster.ToneRaster([30:36],:));
    
    zeroPeakVal = zeroPeak - zeroBase;
    endPeakVal = endPeak-endBase;
    bigStore(1,incInd:incInd + length(s.MBED.HiTrials) - 1) = zeroPeakVal(s.MBED.HiTrials);
    bigStore(2,incInd:incInd + length(s.MBED.LowTrials) - 1) = zeroPeakVal(s.MBED.LowTrials);
    bigStore(3,incInd:incInd + length(s.MBED.HiTrials) - 1) = endPeak(s.MBED.HiTrials);
    bigStore(4,incInd:incInd + length(s.MBED.LowTrials) - 1) = endPeak(s.MBED.LowTrials);
    
    %LETS ALSO PULL LICKS!!!
    lickLatStore = ones(length(s.MBED.ToneDelivery),1);
    for i = 1:length(s.MBED.ToneDelivery)
        lickStore(i) = length(find(s.Licking.ToneRaster(:,2) == i & s.Licking.ToneRaster(:,1) >=0 & s.Licking.ToneRaster(:,1) < 1));
        if find(s.Licking.ToneRaster(:,2) == i & s.Licking.ToneRaster(:,1) >=0,1,'first');
            lickLatStore(i) = s.Licking.ToneRaster(find(s.Licking.ToneRaster(:,2) == i & s.Licking.ToneRaster(:,1) >=0,1,'first'),1);
        else
            lickLatStore(i) = NaN;
        end
        
    end
    bigStore(5,incInd:incInd + length(s.MBED.LowTrials) - 1) = lickStore(s.MBED.HiTrials);
    bigStore(6,incInd:incInd + length(s.MBED.LowTrials) - 1) = lickStore(s.MBED.LowTrials);
    %note that these are with 100ms bins?
    bigStore(7,incInd:incInd + length(s.MBED.LowTrials) - 1) = lickLatStore(s.MBED.HiTrials);
    bigStore(8,incInd:incInd + length(s.MBED.LowTrials) - 1) = lickLatStore(s.MBED.LowTrials);
    %now store velocity data too
    bigStore(9,incInd:incInd + length(s.MBED.LowTrials) - 1) = preVel(s.MBED.HiTrials);
    bigStore(10,incInd:incInd + length(s.MBED.LowTrials) - 1) = preVel(s.MBED.LowTrials);
    bigStore(11,incInd:incInd + length(s.MBED.LowTrials) - 1) = zeroVel(s.MBED.HiTrials);
    bigStore(12,incInd:incInd + length(s.MBED.LowTrials) - 1) = zeroVel(s.MBED.LowTrials);
    
    bigVelStore{bigInd} = s.Locomotion.Velocity;
    photoStore{bigInd} = [s.Photo.Photo.x70dFTime',(s.Photo.Photo.x70dF - preMean)/preSTD];
    
    %now lets also store average photometry traces
    bigPhotAverage(:,bigInd) = mean(zRaster(:,s.MBED.HiTrials)');
    smallPhotAverage(:,bigInd) = mean(zRaster(:,s.MBED.LowTrials)');
    
    %pull licking histograms
    if max(s.Licking.ToneHistHi) > 30
        lickHistHi(:,bigInd) = s.Licking.ToneHistHi/10;
        lickHistLow(:,bigInd) = s.Licking.ToneHistLow/10;
    else
        lickHistHi(:,bigInd) = s.Licking.ToneHistHi;
        lickHistLow(:,bigInd) = s.Licking.ToneHistLow;
    end
%     lickHistLow(31,bigInd) = 0;
    %pull velocities
    velHistHi(:,bigInd) = mean(s.VelRaster.ToneRaster(:,s.MBED.HiTrials)');
    velHistLow(:,bigInd) = mean(s.VelRaster.ToneRaster(:,s.MBED.LowTrials)');
    
    %calculate locomotion ROC over whole trace
    
    smoothVel = smooth(bigVelStore{1}(:,2),5);
    velTrueTime = interp1(s.MBED.Jitter/1000,s.Photo.Jitter,bigVelStore{1}(:,1));
    minVelTime = min(velTrueTime);
    maxVelTime = max(velTrueTime);
    newVelVector = [velTrueTime,smoothVel];
    %remove nan values
    nanFind = find(isnan(newVelVector(:,1)));
    newVelVector(nanFind,:) = [];
    
    peakTimes = s.Photo.Peaks(:,4);
    peakTimes(peakTimes< minVelTime) = [];
    peakTimes(peakTimes > maxVelTime) = [];
    
    [funcOut] = functionLocomotionROC(peakTimes,newVelVector);
    peakAUCStore{bigInd} = funcOut;
    
    %now lets interpolate the photometry signal so they sample at the same rate
    %as velocity
    interpPhot = interp1(photoStore{1}(:,1),photoStore{1}(:,2),newVelVector(:,1));

    minPhot = min(interpPhot);
    maxPhot = max(interpPhot);
    rateInc = 20;

    locomotionInd = double(newVelVector(:,2)>1);
    trueLocs = length(find(locomotionInd == 1));
    truePause = length(find(locomotionInd == 0));

    %rates at which I will threshold as classifier
    rateRange = minPhot:(maxPhot-minPhot)/rateInc:maxPhot;

    rocStore = zeros(4,rateInc);

    for i = 1:rateInc
        %need to conver i to the targeted rate
        threshRate = rateRange(i);
        %find all points at which we classify as locomotion
        classifyInd = (double(interpPhot>=threshRate)+1)*2;
        %compare by adding to locomotionInd
        testInd = classifyInd + locomotionInd;
        %find points with locomotion and classifier(4+1 = 5)
        rocStore(i,1) = length(find(testInd == 5));
        %find points with locomotion but no classifier (2+1 = 3)
        rocStore(i,2) = length(find(testInd == 3));
        %find points with no locomotion but classifier (4 + 0 = 4)
        rocStore(i,3) = length(find(testInd == 4));
        %find points with no locomotion and no classifier (2 + 0 = 2)
        rocStore(i,4) = length(find(testInd == 2));
    end
    %calculate AUC using trapz
    falsePos = rocStore(:,3)/truePause;
    truePos = rocStore(:,1)/trueLocs;
    %reorder in order of false positive from 0 to 1
    [B,I] = sort(falsePos);
    falsePos = B;
    truePos = truePos(I);
    %eliminate duplicate values. 
    [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
    truePos = truePos(ia);
    falsePos = C;
    %calculate estimate of area under curve. 
    AUCoverall(bigInd)= trapz(falsePos,truePos);
    %now lets calculate shuffled AUCs
    numShuff = 1000;
    for shuffInd = 1:numShuff
        shuffleTimes = randperm(length(interpPhot));
        shuffleTimes = interpPhot(shuffleTimes);
        minPhot = min(shuffleTimes);
        maxPhot = max(shuffleTimes);
        rateRange = minPhot:(maxPhot-minPhot)/rateInc:maxPhot;
        rocStore = zeros(4,rateInc);

        for i = 1:rateInc
            %need to conver i to the targeted rate
            threshRate = rateRange(i);
            %find all points at which we classify as locomotion
            classifyInd = (double(shuffleTimes>=threshRate)+1)*2;
            %compare by adding to locomotionInd
            testInd = classifyInd + locomotionInd;
            %find points with locomotion and classifier(4+1 = 5)
            rocStore(i,1) = length(find(testInd == 5));
            %find points with locomotion but no classifier (2+1 = 3)
            rocStore(i,2) = length(find(testInd == 3));
            %find points with no locomotion but classifier (4 + 0 = 4)
            rocStore(i,3) = length(find(testInd == 4));
            %find points with no locomotion and no classifier (2 + 0 = 2)
            rocStore(i,4) = length(find(testInd == 2));
        end
        %calculate AUC using trapz
        falsePos = rocStore(:,3)/truePause;
        truePos = rocStore(:,1)/trueLocs;
        %reorder in order of false positive from 0 to 1
        [B,I] = sort(falsePos);
        falsePos = B;
        truePos = truePos(I);
        %eliminate duplicate values. 
        [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
        truePos = truePos(ia);
        falsePos = C;
        %calculate estimate of area under curve. 
        aucShuffLoco(shuffInd,bigInd)= trapz(falsePos,truePos);
        
        
    end
    
    %now lets do ROC for big vs small trials in terms of licking
    minLick = min(lickStore);
    maxLick = max(lickStore);
    rateRange = minLick:(maxLick-minLick)/rateInc:maxLick;

    rocStore = zeros(4,rateInc);
    
    lickTrues = NaN(1,length(lickStore));
    lickTrues(s.MBED.HiTrials) = 1;
    lickTrues(s.MBED.LowTrials) = 0;

    for i = 1:rateInc
        %need to conver i to the targeted rate
        threshRate = rateRange(i);
        %find all points at which we classify as locomotion
        classifyInd = (double(lickStore>=threshRate)+1)*2;
        %compare by adding to locomotionInd
        testInd = classifyInd + lickTrues;
        %find points with locomotion and classifier(4+1 = 5)
        rocStore(i,1) = length(find(testInd == 5));
        %find points with locomotion but no classifier (2+1 = 3)
        rocStore(i,2) = length(find(testInd == 3));
        %find points with no locomotion but classifier (4 + 0 = 4)
        rocStore(i,3) = length(find(testInd == 4));
        %find points with no locomotion and no classifier (2 + 0 = 2)
        rocStore(i,4) = length(find(testInd == 2));
    end
    
    trueHi = length(s.MBED.HiTrials);
    trueLow = length(s.MBED.LowTrials);

    %calculate AUC using trapz
    falsePos = rocStore(:,3)/trueLow;
    truePos = rocStore(:,1)/trueHi;
    %reorder in order of false positive from 0 to 1
    [B,I] = sort(falsePos);
    falsePos = B;
    truePos = truePos(I);
    %eliminate duplicate values. 
    [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
    truePos = truePos(ia);
    falsePos = C;
    %calculate estimate of area under curve. 
    lickAUC(bigInd)= trapz(falsePos,truePos);
    %do bootstrapping/shuffling
    numShuff = 1000;
    for shuffInd = 1:numShuff
        shuffleTimes = randperm(length(lickStore));
        shuffleTimes = lickStore(shuffleTimes);
        minPhot = min(shuffleTimes);
        maxPhot = max(shuffleTimes);
        rateRange = minPhot:(maxPhot-minPhot)/rateInc:maxPhot;
        rocStore = zeros(4,rateInc);

        for i = 1:rateInc
            %need to conver i to the targeted rate
            threshRate = rateRange(i);
            %find all points at which we classify as locomotion
            classifyInd = (double(shuffleTimes>=threshRate)+1)*2;
            %compare by adding to locomotionInd
            testInd = classifyInd + lickTrues;
            %find points with locomotion and classifier(4+1 = 5)
            rocStore(i,1) = length(find(testInd == 5));
            %find points with locomotion but no classifier (2+1 = 3)
            rocStore(i,2) = length(find(testInd == 3));
            %find points with no locomotion but classifier (4 + 0 = 4)
            rocStore(i,3) = length(find(testInd == 4));
            %find points with no locomotion and no classifier (2 + 0 = 2)
            rocStore(i,4) = length(find(testInd == 2));
        end
        %calculate AUC using trapz
        falsePos = rocStore(:,3)/trueLow;
        truePos = rocStore(:,1)/trueHi;
        %reorder in order of false positive from 0 to 1
        [B,I] = sort(falsePos);
        falsePos = B;
        truePos = truePos(I);
        %eliminate duplicate values. 
        [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
        truePos = truePos(ia);
        falsePos = C;
        %calculate estimate of area under curve. 
        aucShuffLick(shuffInd,bigInd)= trapz(falsePos,truePos);
    end
    %increment index for low trials
    incInd = incInd + length(s.MBED.LowTrials);
    
    
    
    
end

store.bigPhotAverage = bigPhotAverage;
store.smallPhotAverage = smallPhotAverage;
%lets try breaking things into 20 trial blocks

forReps = length(bigStore)/20;
for i = 1:forReps
    condensedBig(:,i) = nanmean(bigStore(:,(i-1)*20+1:(i)*20),2);
    condensedBigVar(:,i) = std(bigStore(:,(i-1)*20+1:(i)*20),0,2);
    condensedSig(1,i) = ranksum(bigStore(1,(i-1)*20+1:(i)*20),bigStore(2,(i-1)*20+1:(i)*20));
    condensedSig(2,i) = ranksum(bigStore(3,(i-1)*20+1:(i)*20),bigStore(4,(i-1)*20+1:(i)*20));
    condensedSig(3,i) = ranksum(bigStore(5,(i-1)*20+1:(i)*20),bigStore(6,(i-1)*20+1:(i)*20));
%     condensedSig(
end

for i = 1:forReps/5
    dayBig(:,i) = nanmean(bigStore(:,(i-1)*100+1:(i)*100),2);
    dayVar(:,i) = std(bigStore(:,(i-1)*100+1:(i)*100),0,2);
    daySig(1,i) = ranksum(bigStore(1,(i-1)*100+1:(i)*100),bigStore(2,(i-1)*100+1:(i)*100));
    daySig(2,i) = ranksum(bigStore(3,(i-1)*100+1:(i)*100),bigStore(4,(i-1)*100+1:(i)*100));
    daySig(3,i) = ranksum(bigStore(5,(i-1)*100+1:(i)*100),bigStore(6,(i-1)*100+1:(i)*100));
end

bigStore(1,bigStore(1,:) <0) = 0;
bigStore(2,bigStore(2,:) <0) = 0;
bigStore(9,:) = (bigStore(1,:)-bigStore(2,:))./(bigStore(1,:)+bigStore(2,:));
%fix NANs from 0/0
bigStore(9,isnan(bigStore(9,:))) = 0;

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);


zeroPeakVel = zeroPeak - zeroBase;

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
plot(bigStore(1,:),bigStore(11,:)-bigStore(9,:),'b.')
hold on
plot(bigStore(2,:),bigStore(12,:)-bigStore(10,:),'r.')
title('Scatter of Photometry Peak with Mean Velocity')
xlabel('Photometry Peak Size')
ylabel('Vel Change Relative to Baseline')
savefig(hFig,'vel and photScatter');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'vel and photScatter','-dpdf','-r0')



%now also plot out average photometry traces
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

photoAxis = [-3:0.021:5.001];
photMax = max([max(max(bigPhotAverage)),max(max(smallPhotAverage))]);
photMin = min([min(min(bigPhotAverage)),min(min(smallPhotAverage))]);

lickMax = max([max(max(lickHistHi)),max(max(1))]);
lickMin = min([min(min(lickHistHi)),min(min(1))]);

velMax = max([max(max(velHistHi)),max(max(velHistLow))]);
velMin = min([min(min(velHistHi)),min(min(velHistLow))]);

subplot(3,2,1) %plot DS
hold on
for i = 1:numFiles
    plot(photoAxis,bigPhotAverage(:,i),'LineWidth',2,'Color',[i/numFiles 0 0])
    if i == 1 | i == numFiles
        plot(photoAxis,bigPhotAverage(:,i),'LineWidth',4,'Color',[i/numFiles 0 0])
    end
end
ylim([photMin photMax])
xlim([-3 5])
title('DS Photometry')

subplot(3,2,2) %plot NS
hold on
for i = 1:numFiles
    plot(photoAxis,smallPhotAverage(:,i),'LineWidth',2,'Color',[i/numFiles 0 0])
    if i == 1 | i == numFiles
        plot(photoAxis,smallPhotAverage(:,i),'LineWidth',4,'Color',[i/numFiles 0 0])
    end
end
ylim([photMin photMax])
xlim([-3 5])
title('NS Photometry')

subplot(3,2,3) %plot DS velocity trace
hold on
for i = 1:numFiles
    plot([-3:0.1:5],velHistHi(:,i),'LineWidth',2,'Color',[i/numFiles 0 0])
    if i == 1 | i == numFiles
        plot([-3:0.1:5],velHistHi(:,i),'LineWidth',4,'Color',[i/numFiles 0 0])
    end
end
ylim([velMin velMax])
xlim([-3 5])
title('DS VELOCITY')

subplot(3,2,4) %plot NS velocity trace
hold on
for i = 1:numFiles
    plot([-3:0.1:5],velHistLow(:,i),'LineWidth',2,'Color',[i/numFiles 0 0])
    if i == 1 | i == numFiles
        plot([-3:0.1:5],velHistLow(:,i),'LineWidth',4,'Color',[i/numFiles 0 0])
    end
end
ylim([velMin velMax])
xlim([-3 5])
title('NS VELOCITY')

subplot(3,2,5) %plot DS Licking trace
hold on
for i = 1:numFiles
    plot([-3:0.1:5],lickHistHi(:,i),'LineWidth',2,'Color',[i/numFiles 0 0])
    if i == 1 | i == numFiles
        plot([-3:0.1:5],lickHistHi(:,i),'LineWidth',4,'Color',[i/numFiles 0 0])
    end
end
ylim([lickMin lickMax])
xlim([-3 5])
title('DS Licking')

subplot(3,2,6) %plot NS licking trace
hold on
for i = 1:numFiles
    plot([-3:0.1:5],lickHistLow(:,i),'LineWidth',2,'Color',[i/numFiles 0 0])
    if i == 1 | i == numFiles
        plot([-3:0.1:5],lickHistLow(:,i),'LineWidth',4,'Color',[i/numFiles 0 0])
    end
end
ylim([lickMin lickMax])
xlim([-3 5])
title('NS Licking')

savefig(hFig,'photVSlickVSvel');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'photVSlickVSvel','-dpdf','-r0')


hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
hold on
plot(bigStore(1,:),'bo')
plot(smooth(bigStore(1,:),11,'lowess'),'b.-')
plot(bigStore(2,:),'ro')
plot(smooth(bigStore(2,:),11,'lowess'),'r.-')
for i = 1:(numFiles)
    plot([i*100 i*100],[0 0.1],'k')
end
xlim([0 length(bigStore)])
title('DS(b) and NS(r) ONSET Plotted by Trial, 11 trial moving average.')

subplot(2,2,3)
hold on
plot(bigStore(9,:),'ko')
plot(smooth(bigStore(9,:),11,'lowess'),'k.-')
plot([1 length(bigStore)],[1 1],'g','LineWidth',2)
for i = 1:(numFiles)
    plot([(i-1)*100+1 i*100],[mean(bigStore(9,[(i-1)*100+1:i*100])) mean(bigStore(9,[(i-1)*100+1:i*100]))],'r','LineWidth',2)
    plot([i*100 i*100],[-1 1],'k')
end

xlim([0 length(bigStore)])
title('DS(b)/NS(r) Plotted by Trial, 11 trial moving average.')

%licks!
subplot(2,2,2)
hold on
plot(bigStore(5,:),'bo')
plot(smooth(bigStore(5,:),11,'lowess'),'b.-')
plot(bigStore(6,:),'ro')
plot(smooth(bigStore(6,:),11,'lowess'),'r.-')

for i = 1:(numFiles)
    plot([i*100 i*100],[0 7],'k')
end
xlim([0 length(bigStore)])
title('DS(b) and NS(r) Licking Plotted by Trial, 11 trial moving average.')


%PLOT OUT LATENCY
subplot(2,2,4)
hold on
plot(bigStore(7,:),'bo')
plot(smooth(bigStore(7,:),11,'lowess'),'b.-')
plot(bigStore(8,:),'ro')
plot(smooth(bigStore(8,:),11,'lowess'),'r.-')
for i = 1:(numFiles)
    plot([i*100 i*100],[0 7],'k')
end

xlim([0 length(bigStore)])
ylim([0 2])
title('DS(b) and NS(r)  Licking Plotted by Trial, 11 trial moving average.')

savefig(hFig,'dFoF vs Licking');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'dFoF vs Licking','-dpdf','-r0')








hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(4,2,1)
hold on
plot(condensedBig(1,:),'b.-')
plot(condensedBig(2,:),'r.-')
findSig = find(condensedSig(1,:) < 0.05);
plot(findSig,condensedBig(1,findSig),'g*')
for i = 1:(numFiles)
    plot([i*5 i*5],[0 0.1],'k')
end
xlim([0 length(condensedBig)])
title('DS(b) and NS(r) ONSET Plotted by Trial,binned 20trial average')

subplot(4,2,3)
hold on
plot(condensedBigVar(1,:),'b.-')
plot(condensedBigVar(2,:),'r.-')
for i = 1:(numFiles)
    plot([i*5 i*5],[0 0.1],'k')
end
xlim([0 length(condensedBig)])
title('DS(b) and NS(r) VAR Plotted by Trial,binned 20trial average')

subplot(4,2,5)
hold on
plot(condensedBig(5,:),'b.-')
plot(condensedBig(6,:),'r.-')
for i = 1:(numFiles)
    plot([i*5 i*5],[0 7],'k')
end
xlim([0 length(condensedBig)])
title('DS(b) and NS(r) Licking Plotted by Trial, binned 20trial average')

prefScore = (condensedBig(5,:)) ./ (condensedBig(5,:) + condensedBig(6,:));

subplot(4,2,7)
hold on
plot(prefScore,'g.-')
findSig = find(condensedSig(3,:) < 0.05);
plot(findSig,prefScore(findSig),'ro')
plot([0 length(condensedBig)],[0 0],'k')
for i = 1:(numFiles)
    plot([i*5 i*5],[-1 1],'k')
end

plot([0 length(condensedBig)],[behavLim behavLim],'r')
ylim([0 1])
xlim([0 length(condensedBig)])
title('Licking Preference Plotted by Trial, binned 20trial average')

subplot(4,2,2)
hold on
plot(dayBig(1,:),'b.-')
plot(dayBig(2,:),'r.-')
findSig = find(daySig(1,:) < 0.05);
plot(findSig,dayBig(1,findSig),'g*')

xlim([0 length(dayBig)])
title('DS(b) and NS(r) ONSET Plotted by DAY')

subplot(4,2,4)
hold on
plot(dayVar(1,:),'b.-')
plot(dayVar(2,:),'r.-')

xlim([0 length(dayBig)])
title('DS(b) and NS(r) VAR Plotted by DAY')

subplot(4,2,6)
hold on
plot(dayBig(5,:),'b.-')
plot(dayBig(6,:),'r.-')

xlim([0 length(dayBig)])
title('DS(b) and NS(r) Licking Plotted by DAY')

% prefScore = (dayBig(5,:) - dayBig(6,:)) ./ (dayBig(5,:) + dayBig(6,:));
prefScore = (dayBig(5,:)) ./ (dayBig(5,:) + dayBig(6,:));

subplot(4,2,8)
hold on
plot(prefScore,'g.-')
findSig = find(daySig(3,:) < 0.05);
plot(findSig,prefScore(findSig),'ro')
plot([0 length(dayBig)],[0 0],'k')
ylim([0 1])
plot([0 length(dayBig)],[behavLim behavLim],'r')
xlim([0 length(dayBig)])
title('Licking Preference Plotted by DAY')

savefig(hFig,'dFoF vs Licking');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'condensed dFoF vs Licking','-dpdf','-r0')


%now lets extract data from the relevant time points. I want to get all the
%info at several specific times: day 1, which should be early, the first
%day of at least two consecutive days in which the animal performs above
%criterion, and the last day of behavior

%to find the first day above criterion, lets look at hte day to day data. 
whileTrig = 0;
setTrig = 0;
indStart = 1;
while whileTrig == 0
    %check at the current index
    prefVal = prefScore(indStart);
    if prefVal >= behavLim & setTrig == 0
        setTrig = 1;
        disp('First Threshold Crossing')
        indStart = indStart + 1;
    elseif prefVal >= behavLim & setTrig == 1
        disp('TargetFound')
        targetInd = indStart -1;
        whileTrig = 1;
        disp(strcat('DAY ',num2str(targetInd)))
        break
    elseif prefVal < behavLim & setTrig == 1
        disp('Failure To Maintain Behavior')
        setTrig = 0;
        indStart = indStart + 1;
    elseif prefVal < behavLim & setTrig == 0
        disp('No Threshold Crossing')
        indStart = indStart + 1;
    end
        
end




%now lets try calculation ROC (again Q_Q)

rateInc = 20;
locomotionInd = repmat([1 0],1,100);
trueLocs = length(find(locomotionInd == 1));
truePause = length(find(locomotionInd == 0));

%find peaks from average histograms
for ind = 1:numFiles
    smoothRate = reshape(bigStore([1:2],[1+100*(ind-1):100*ind]),1,[]);
    %find min and max rates!
    minRate = (min(smoothRate));
    maxRate = (max(smoothRate));

    %rates at which I will threshold as classifier
    rateRange = minRate:(maxRate-minRate)/rateInc:maxRate;
    
    rocStore = zeros(4,rateInc);

    for i = 1:rateInc
        %need to conver i to the targeted rate
        threshRate = rateRange(i);
        %find all points at which we classify as locomotion
        classifyInd = (double(smoothRate>=threshRate)+1)*2;
        %compare by adding to locomotionInd
        testInd = classifyInd + locomotionInd;
        %find points with locomotion and classifier(4+1 = 5)
        rocStore(i,1) = length(find(testInd == 5));
        %find points with locomotion but no classifier (2+1 = 3)
        rocStore(i,2) = length(find(testInd == 3));
        %find points with no locomotion but classifier (4 + 0 = 4)
        rocStore(i,3) = length(find(testInd == 4));
        %find points with no locomotion and no classifier (2 + 0 = 2)
        rocStore(i,4) = length(find(testInd == 2));
    end
    
    %calculate AUC using trapz
    falsePos = rocStore(:,3)/truePause;
    truePos = rocStore(:,1)/trueLocs;
    %reorder in order of false positive from 0 to 1
    [B,I] = sort(falsePos);
    falsePos = B;
    truePos = truePos(I);
    %eliminate duplicate values. 
    [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
    truePos = truePos(ia);
    falsePos = C;
    %calculate estimate of area under curve. 
    trueAUC(ind)= trapz(falsePos,truePos);
end

%now calculate AUC prediction of velocity relative to photometry



%now look at AUC for locomotion predicting trial type
locoChange = reshape(bigStore([11:12],:) - bigStore([9:10],:),1,[]);

rateInc = 20;
locomotionInd = repmat([1 0],1,100);
trueLocs = length(find(locomotionInd == 1));
truePause = length(find(locomotionInd == 0));
try
    %find peaks from average histograms
    for ind = 1:numFiles
        smoothRate = locoChange(:,[1+200*(ind-1):200*ind]);
        %find min and max rates!
        minRate = (min(smoothRate));
        maxRate = (max(smoothRate));

        %rates at which I will threshold as classifier
        rateRange = minRate:(maxRate-minRate)/rateInc:maxRate;

        rocStore = zeros(4,rateInc);

        for i = 1:rateInc
            %need to conver i to the targeted rate
            threshRate = rateRange(i);
            %find all points at which we classify as locomotion
            classifyInd = (double(smoothRate>=threshRate)+1)*2;
            %compare by adding to locomotionInd
            testInd = classifyInd + locomotionInd;
            %find points with locomotion and classifier(4+1 = 5)
            rocStore(i,1) = length(find(testInd == 5));
            %find points with locomotion but no classifier (2+1 = 3)
            rocStore(i,2) = length(find(testInd == 3));
            %find points with no locomotion but classifier (4 + 0 = 4)
            rocStore(i,3) = length(find(testInd == 4));
            %find points with no locomotion and no classifier (2 + 0 = 2)
            rocStore(i,4) = length(find(testInd == 2));
        end

        %calculate AUC using trapz
        falsePos = rocStore(:,3)/truePause;
        truePos = rocStore(:,1)/trueLocs;
        %reorder in order of false positive from 0 to 1
        [B,I] = sort(falsePos);
        falsePos = B;
        truePos = truePos(I);
        %eliminate duplicate values. 
        [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
        truePos = truePos(ia);
        falsePos = C;
        %calculate estimate of area under curve. 
        locoChangeAUC(ind)= trapz(falsePos,truePos);
    end
catch
    locoChangeAUC = zeros(numFiles,1);
    disp('Loco Change AUC Analysis Failed')
end

try
    %just use mean values at zero point, try again
    locoChange = reshape(bigStore([11:12],:),1,[]);

    rateInc = 20;
    locomotionInd = repmat([1 0],1,100);
    trueLocs = length(find(locomotionInd == 1));
    truePause = length(find(locomotionInd == 0));

    %find peaks from average histograms
    for ind = 1:numFiles
        smoothRate = locoChange(:,[1+200*(ind-1):200*ind]);
        %find min and max rates!
        minRate = (min(smoothRate));
        maxRate = (max(smoothRate));

        %rates at which I will threshold as classifier
        rateRange = minRate:(maxRate-minRate)/rateInc:maxRate;

        rocStore = zeros(4,rateInc);

        for i = 1:rateInc
            %need to conver i to the targeted rate
            threshRate = rateRange(i);
            %find all points at which we classify as locomotion
            classifyInd = (double(smoothRate>=threshRate)+1)*2;
            %compare by adding to locomotionInd
            testInd = classifyInd + locomotionInd;
            %find points with locomotion and classifier(4+1 = 5)
            rocStore(i,1) = length(find(testInd == 5));
            %find points with locomotion but no classifier (2+1 = 3)
            rocStore(i,2) = length(find(testInd == 3));
            %find points with no locomotion but classifier (4 + 0 = 4)
            rocStore(i,3) = length(find(testInd == 4));
            %find points with no locomotion and no classifier (2 + 0 = 2)
            rocStore(i,4) = length(find(testInd == 2));
        end

        %calculate AUC using trapz
        falsePos = rocStore(:,3)/truePause;
        truePos = rocStore(:,1)/trueLocs;
        %reorder in order of false positive from 0 to 1
        [B,I] = sort(falsePos);
        falsePos = B;
        truePos = truePos(I);
        %eliminate duplicate values. 
        [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
        truePos = truePos(ia);
        falsePos = C;
        %calculate estimate of area under curve. 
        locoValAUC(ind)= trapz(falsePos,truePos);
    end
catch
    disp('Loco Mean Analysis Failed')
    locoValAUC = zeros(numFiles,1);
end
%looks like the actual mean velocity value is more helpful as a predictor
%than the change in velocity?

%LOOK AT LICKING AS AUC
locoChange = reshape(bigStore([5:6],:),1,[]);

rateInc = 20;
locomotionInd = repmat([1 0],1,100);
trueLocs = length(find(locomotionInd == 1));
truePause = length(find(locomotionInd == 0));

%find peaks from average histograms
try
    for ind = 1:numFiles
        smoothRate = locoChange(:,[1+200*(ind-1):200*ind]);
        %find min and max rates!
        minRate = (min(smoothRate));
        maxRate = (max(smoothRate));

        %rates at which I will threshold as classifier
        rateRange = minRate:(maxRate-minRate)/rateInc:maxRate;

        rocStore = zeros(4,rateInc);

        for i = 1:rateInc
            %need to conver i to the targeted rate
            threshRate = rateRange(i);
            %find all points at which we classify as locomotion
            classifyInd = (double(smoothRate>=threshRate)+1)*2;
            %compare by adding to locomotionInd
            testInd = classifyInd + locomotionInd;
            %find points with locomotion and classifier(4+1 = 5)
            rocStore(i,1) = length(find(testInd == 5));
            %find points with locomotion but no classifier (2+1 = 3)
            rocStore(i,2) = length(find(testInd == 3));
            %find points with no locomotion but classifier (4 + 0 = 4)
            rocStore(i,3) = length(find(testInd == 4));
            %find points with no locomotion and no classifier (2 + 0 = 2)
            rocStore(i,4) = length(find(testInd == 2));
        end

        %calculate AUC using trapz
        falsePos = rocStore(:,3)/truePause;
        truePos = rocStore(:,1)/trueLocs;
        %reorder in order of false positive from 0 to 1
        [B,I] = sort(falsePos);
        falsePos = B;
        truePos = truePos(I);
        %eliminate duplicate values. 
        [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
        truePos = truePos(ia);
        falsePos = C;
        %calculate estimate of area under curve. 
        lickAUC(ind)= trapz(falsePos,truePos);
    end
catch
    lickAUC = zeros(numFiles,1);
    disp('LICK AUC ANALYSIS FAILED')
end

hFig = figure;
plot(trueAUC,'k','LineWidth',2)
hold on
plot(AUCoverall,'g')
for i = 1:numFiles
    peakAUCVal(i) = peakAUCStore{i}.TrueAUC;
end
plot(peakAUCVal,'g*-')
plot(lickAUC,'b')
plot(locoChangeAUC,'r')
plot(locoValAUC,'r*-')
legend
title('AUC Calculation for Photometry Kphoto Gloco')
ylabel('AUC Score')
xlabel('Days')


auc.trueAUC = trueAUC;
auc.overallAUC = AUCoverall;
auc.peaksOverall = peakAUCVal;
auc.lickAUC = lickAUC;
auc.locoChangeAUC = locoChangeAUC;
auc.locoValAUC = locoValAUC;


%now lets try for licks

keyVals = zeros(3,1);

keyVals(1,1) = prefScore(1);
keyVals(1,2) = dayBig(1,1);
keyVals(1,3) = dayBig(2,1);
keyVals(1,4) = dayBig(5,1);
keyVals(1,5) = dayBig(6,1);
keyVals(1,6) = AUCoverall(1);
keyVals(1,7) = trueAUC(1);

keyVals(2,1) = prefScore(targetInd);
keyVals(2,2) = dayBig(1,targetInd);
keyVals(2,3) = dayBig(2,targetInd);
keyVals(2,4) = dayBig(5,targetInd);
keyVals(2,5) = dayBig(6,targetInd);
keyVals(2,6) = AUCoverall(targetInd);
keyVals(2,7) = trueAUC(targetInd);

keyVals(3,1) = prefScore(end);
keyVals(3,2) = dayBig(1,end);
keyVals(3,3) = dayBig(2,end);
keyVals(3,4) = dayBig(5,end);
keyVals(3,5) = dayBig(6,end);
keyVals(3,6) = AUCoverall(end);
keyVals(3,7) = trueAUC(end);

figure
hold on
% plot(keyVals(:,1),'g')
plot(keyVals(:,2))
plot(keyVals(:,3),'r')





%try reshaping data to days
% dailyZeroBig = reshape(dayBig(1,:),5,[]);

%hmmm not so useful...dont really see strong patterns there.

%lets compute preference scores with pre-empt licks for correct over total

% prefScore = condensedBig(5,:) ./ (condensedBig(5,:) + condensedBig(6,:));

newFileName = strcat(targetFiles{1}(1:9),'.mat');

save(newFileName,'bigStore','prefScore','condensedBig','condensedSig','keyVals','auc','store');
