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
    %baseline period will simply be average of the five bins before zero
    zeroBase = mean(s.PhotoRaster.ToneRaster(zeroPoint - 4:zeroPoint,:));
    endBase =  mean(s.PhotoRaster.ToneRaster(endPoint - 4:endPoint,:));
    %find peak values for the twentyfive bins following the stimulus onset
    %(~500ms)
    zeroPeak = max(s.PhotoRaster.ToneRaster(zeroPoint:zeroPoint+ 25,:));
    endPeak = max(s.PhotoRaster.ToneRaster(zeroPoint:zeroPoint+ 25,:));
    
    zeroPeakVal = zeroPeak - zeroBase;
    endPeakVal = endPeak-endBase;
    bigStore(1,incInd:incInd + length(s.MBED.HiTrials) - 1) = zeroPeakVal(s.MBED.HiTrials);
    bigStore(2,incInd:incInd + length(s.MBED.LowTrials) - 1) = zeroPeakVal(s.MBED.LowTrials);
    bigStore(3,incInd:incInd + length(s.MBED.HiTrials) - 1) = endPeak(s.MBED.HiTrials);
    bigStore(4,incInd:incInd + length(s.MBED.LowTrials) - 1) = endPeak(s.MBED.LowTrials);
%     bigStore(5,incInd:incInd + length(s.MBED.LowTrials) - 1) = 0;
%     bigStore(5,s.MBED.HiTrials + incInd - 1) = 1;
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
    
    peakValStore{bigInd} = s.Photo.Peaks(:,1);
    
    incInd = incInd + length(s.MBED.LowTrials);
    
    
    
end


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

prefScore = (condensedBig(5,:) - condensedBig(6,:)) ./ (condensedBig(5,:) + condensedBig(6,:));

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

keyVals = zeros(3,1);

keyVals(1,1) = prefScore(1);
keyVals(1,2) = dayBig(1,1);
keyVals(1,3) = dayBig(2,1);
keyVals(1,4) = dayBig(5,1);
keyVals(1,5) = dayBig(6,1);

keyVals(2,1) = prefScore(targetInd);
keyVals(2,2) = dayBig(1,targetInd);
keyVals(2,3) = dayBig(2,targetInd);
keyVals(2,4) = dayBig(5,targetInd);
keyVals(2,5) = dayBig(6,targetInd);

keyVals(3,1) = prefScore(end);
keyVals(3,2) = dayBig(1,end);
keyVals(3,3) = dayBig(2,end);
keyVals(3,4) = dayBig(5,end);
keyVals(3,5) = dayBig(6,end);

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

save(newFileName,'bigStore','prefScore','condensedBig','condensedSig','keyVals');
