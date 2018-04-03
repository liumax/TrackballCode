

fileNames = what;
fileNames = fileNames.mat;
masterIndex = strfind(fileNames,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
fileNames = fileNames(masterIndex);

fullNumFiles = length(fileNames);

fullStore = cell(fullNumFiles,5);


for fullInd = 1:fullNumFiles
    load(fileNames{fullInd})
    %store file name
    fullStore{fullInd,1} = fileNames{fullInd};
    %store preference data
    fullStore{fullInd,2} = prefScore;
    %store AUC data for licking predicting trials
    fullStore{fullInd,3} = lickAUC;
    fullStore{fullInd,4} = findSig;
    %store AUC data for photometry predicting trial type
    fullStore{fullInd,5} = trueAUC;
   
    
    %store information from photometry data
    fullStore{fullInd,7} = bigPhotAverage;
    fullStore{fullInd,8} = smallPhotAverage;
    fullStore{fullInd,9} = velHistHi;
    fullStore{fullInd,10} = velHistLow;
    fullStore{fullInd,11} = lickHistHi;
    fullStore{fullInd,12} = lickHistLow;
    fullStore{fullInd,13} = dayBig(1:2,:);
    fullStore{fullInd,14} = daySig(1,:);
    fullStore{fullInd,15} = dayBig(5:6,:);

end

fullNumFiles = length(fileNames);

%lets plot out preferences for licking!! first extract data
prefHolder = NaN(18,20);
for i = 1:fullNumFiles
    lengthFind = length(fullStore{i,2});
    prefHolder(i,1:lengthFind) = fullStore{i,2};
end

%plot out results
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

hFig = figure;
hold on
for i = 1:fullNumFiles
    plot(prefHolder(i,:),'LineWidth',1)
    findLast = find(prefHolder(i,:)>0,1,'last');
    plot(findLast,prefHolder(i,findLast),'b*')
    %try and find the first above the cutoff. 
end
ylim([0 1])
xlim([1 14])
xlabel('Days of Training')
ylabel('Preference Score')

%now lets look at how well licking works as a ROC

lickHolder = NaN(18,20);
for i = 1:fullNumFiles
    lengthFind = length(fullStore{i,3});
    lickHolder(1:lengthFind,i) = fullStore{i,3};
end

lickSigHolder = NaN(18,20);
for i = 1:fullNumFiles
    lengthFind = length(fullStore{i,4});
    lickSigHolder(1:lengthFind,i) = fullStore{i,4};
end

%plot out results
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

hFig = figure;
hold on
for i = 1:fullNumFiles
    plot(lickHolder(:,i),'LineWidth',1)
    findLast = find(lickHolder(:,i)>0,1,'last');
%     findReal = find(~isnan(lickSigHolder(:,i)));
    plot(findLast,lickHolder(findLast,i),'b*')
%     plot(lickSigHolder(findReal,i),lickHolder(lickSigHolder(findReal,i),i),'r.')
    %try and find the first above the cutoff. 
    %find mean and plot
    tester = nanmean(lickHolder');
    plot(tester,'LineWidth',2)
end
ylim([0 1])
xlim([1 14])
xlabel('Days of Training')
ylabel('AUC Score for Licking')

%plot with different colors
hFig = figure;
plot(lickHolder)
hold on
tester = nanmean(lickHolder');
plot(tester,'LineWidth',2)
ylim([0 1])
xlim([1 14])
xlabel('Days of Training')
ylabel('AUC Score for Licking')
savefig(hFig,'trainingBehaviorAUCProgress');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'trainingBehaviorAUCProgress','-dpdf','-r0')
% print(hFig,'trainingBehaviorAUCProgress','-djpeg','-r0')


%now lets figure out first day of significant behavior followed by further
%successful behavior
firstSig = NaN(fullNumFiles,1);
tester = diff(lickSigHolder);
for i = 1:fullNumFiles
    findFirst = find(tester(:,i) == 1,1,'first');
    if findFirst
        firstSig(i) = lickSigHolder(find(tester(:,i) == 1,1,'first'),i);
    end
end

%insert data about behavior!
targetFreqs = [8000	13454
6727	4000
8000	16000
9513	16000
5657	11314
11314	6727
19027	11314
6727	11314
9513	5657
8000	13454
13454	8000
8000	16000
11314	6727
13454	9514
8000	11314
9514	13454
11314	8000
8000	11314];


%do best fit line
useValue = find(~isnan(firstSig));
% p = polyfit(spectSep(useValue),firstSig(useValue)',1);
% y1 = polyval(p,targetFreqs(useValue,1));

hFig = figure
semilogx(targetFreqs(useValue,1),firstSig(useValue),'r.')
% hold on
% plot(targetFreqs(useValue,1),y1);
xlim([4000 20000])
xlabel('Frequency (Hz)')
ylabel('Days to Significant AUC')
savefig(hFig,'targetFreqVsLearnSpeed');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'targetFreqVsLearnSpeed','-dpdf','-r0')
%dont see really strong effect of target frequency on learning. maybe
%spectral separation?

for i = 1:fullNumFiles
    smallFreq = min(targetFreqs(i,:));
    bigFreq = max(targetFreqs(i,:));
    spectSep(i) = log(bigFreq/smallFreq)/log(2);
end


%do best fit line
% p = polyfit(spectSep(useValue),firstSig(useValue)',1);
% y1 = polyval(p,spectSep);

hFig= figure
plot(spectSep,firstSig,'r.')
xlabel('Spectral Separation (Octaves)')
ylabel('Days to Significant AUC')
savefig(hFig,'spectSeparationVsLearnSpeed');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'spectSeparationVsLearnSpeed','-dpdf','-r0')
%looks like spectral separation is bigger determinant of behavior

%plot histogram of first significant day
firstmin = min(firstSig);
firstmax = max(firstSig);

hFig = figure
hist(firstSig,[firstmin:firstmax])
xlim([firstmin firstmax])
xlabel('Days to First Significant AUC')
ylabel('Number of Mice')
title('Histogram of Days to Significant Behavioral AUC')
savefig(hFig,'sigHist');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'sigHist','-dpdf','-r0')



%% now what do we do? Lets plot out course of photometry changes

d1 = [1,2,5,6,7,8,9,16,17,18];
a2a = [3,4,10,11,12,13,14,15];


%need to plot out d1 photometry peak values
hFig = figure
hold on
DSstore = NaN(20,length(d1));
NSstore = NaN(20,length(d1));
for i = 1:length(d1)
    lengthFind = length(fullStore{d1(i),13}(1,:));
    DSstore(1:lengthFind,i) = fullStore{d1(i),13}(1,:);
    NSstore(1:lengthFind,i) = fullStore{d1(i),13}(2,:);
    plot(fullStore{d1(i),13}(1,:),'b.-')
    plot(fullStore{d1(i),13}(2,:),'r.-')
end
DSMean = nanmean(DSstore');
NSMean = nanmean(NSstore');
plot(DSMean,'b','LineWidth',2)
plot(NSMean,'r','LineWidth',2)
xlim([1 14])
ylabel('Z-Scored Peak Response')
xlabel('Days of Training')
title('Change in Daily D1 Peak Photometry Response Over Days')
savefig(hFig,'AllD1ZScorePhotometry');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AllD1ZScorePhotometry','-dpdf','-r0')

%now lets try normalizing to the first day

hFig = figure
hold on
normDS = NaN(20,length(d1));
normNS = NaN(20,length(d1));
for i = 1:length(d1)
%     lengthFind = length(fullStore{d1(i),13}(1,:));
    normDS(:,i) = DSstore(:,i)/DSstore(1,i);
    normNS(:,i) = NSstore(:,i)/NSstore(1,i);
    plot(DSstore(:,i)/DSstore(1,i),'b.-')
    plot(NSstore(:,i)/NSstore(1,i),'r.-')
end
DSMean = nanmean(normDS');
NSMean = nanmean(normNS');
plot(DSMean,'b','LineWidth',2)
plot(NSMean,'r','LineWidth',2)
xlim([1 14])
ylim([0 3])
ylabel('Normalized Peak Response')
xlabel('Days of Training')
title('Normalized Change in Daily D1 Peak Photometry Response Over Days')
% %not sure I like this better, keep the other one better!
savefig(hFig,'AllD1ZScorePhotometryNorm');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AllD1ZScorePhotometryNorm','-dpdf','-r0')

%now lets simply plot as single points of first and last session
firstLastStore = zeros(length(d1),4);
for i = 1:length(d1)
    firstLastStore(i,1) = DSstore(1,i);
    firstLastStore(i,3) = NSstore(1,i);
    %find the last numeric value
    findNums = find(~isnan(DSstore(:,i)));
    firstLastStore(i,2) = DSstore(findNums(end),i);
    firstLastStore(i,4) = NSstore(findNums(end),i);
end

hFig = figure;
hold on
for i = 1:length(d1)
    plot(firstLastStore(i,1:2),'b.-')
    plot(firstLastStore(i,3:4),'r.-')
end
errorbar(mean(firstLastStore(:,1:2)),std(firstLastStore(:,1:2))/sqrt(length(d1)),'b.-','LineWidth',2)
errorbar(mean(firstLastStore(:,3:4)),std(firstLastStore(:,3:4))/sqrt(length(d1)),'r.-','LineWidth',2)
xlim([0 3])
xlabel('First vs Last Day')
ylabel('Photometry Signal (Zscored)')
title('D1 First and Last Day')
savefig(hFig,'d1firstVSlastPhotValue');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'d1firstVSlastPhotValue','-dpdf','-r0')

%now we can do ratio of change between the two
firstLastRatio = zeros(length(d1),2);
firstLastRatio(:,1) = firstLastStore(:,2)./firstLastStore(:,1);
firstLastRatio(:,2) = firstLastStore(:,4)./firstLastStore(:,3);

hFig = figure;
hold on
for i = 1:length(d1)
    plot(firstLastRatio(i,:),'k')
end
plot([ones(length(d1),1)],firstLastRatio(:,1),'b.')
plot([ones(length(d1),1)*2],firstLastRatio(:,2),'r.')
xlim([0 3])
xlabel('DS (1) and NS(2)')
ylabel('Last/First Day Photometry')
title('D1 Ratio Change First and Last Day')
savefig(hFig,'d1firstVSlastPhotRatio');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'d1firstVSlastPhotRatio','-dpdf','-r0')

%we can also try linear regression fit for all the sets of behavioral data

for i = 1:length(d1)
    x = find(~isnan(DSstore(:,i)));
    y1 = DSstore(x,i);
    P = polyfit(x,y1,1);
    slopeStore(i,1) = P(1);
    y1 = NSstore(x,i);
    P = polyfit(x,y1,1);
    slopeStore(i,2) = P(1);
end


hFig = figure;
hold on
for i = 1:length(d1)
    plot(slopeStore(i,:),'k')
end
plot([ones(length(d1),1)],slopeStore(:,1),'b.')
plot([ones(length(d1),1)*2],slopeStore(:,2),'r.')
xlim([0 3])
xlabel('DS (1) and NS(2)')
ylabel('Slope for Overall Behavior')
title('D1 Slope First and Last Day')
savefig(hFig,'d1firstVSlastPhotSlope');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'d1firstVSlastPhotSlope','-dpdf','-r0')

%% now do A2A
hFig = figure
hold on
DSstore = NaN(20,length(a2a));
NSstore = NaN(20,length(a2a));
for i = 1:length(a2a)
    lengthFind = length(fullStore{a2a(i),13}(1,:));
    DSstore(1:lengthFind,i) = fullStore{a2a(i),13}(1,:);
    NSstore(1:lengthFind,i) = fullStore{a2a(i),13}(2,:);
    plot(fullStore{a2a(i),13}(1,:),'b.-')
    plot(fullStore{a2a(i),13}(2,:),'r.-')
end
DSMean = nanmean(DSstore');
NSMean = nanmean(NSstore');
plot(DSMean,'b','LineWidth',2)
plot(NSMean,'r','LineWidth',2)
xlim([1 14])
ylabel('Z-Scored Peak Response')
xlabel('Days of Training')
title('Change in Daily A2A Peak Photometry Response Over Days')
savefig(hFig,'AllA2AZScorePhotometry');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AllA2AZScorePhotometry','-dpdf','-r0')

hFig = figure
hold on
normDS = NaN(20,length(a2a));
normNS = NaN(20,length(a2a));
for i = 1:length(a2a)
%     lengthFind = length(fullStore{d1(i),13}(1,:));
    normDS(:,i) = DSstore(:,i)/DSstore(1,i);
    normNS(:,i) = NSstore(:,i)/NSstore(1,i);
    plot(DSstore(:,i)/DSstore(1,i),'b.-')
    plot(NSstore(:,i)/NSstore(1,i),'r.-')
end
DSMean = nanmean(normDS');
NSMean = nanmean(normNS');
plot(DSMean,'b','LineWidth',2)
plot(NSMean,'r','LineWidth',2)
xlim([1 14])
ylabel('Normalized Peak Response')
xlabel('Days of Training')
title('Normalized Change in Daily A2A Peak Photometry Response Over Days')
% %not sure I like this better, keep the other one better!
savefig(hFig,'AllA2AZScorePhotometryNorm');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'AllA2AZScorePhotometryNorm','-dpdf','-r0')


%now lets simply plot as single points of first and last session
firstLastStore = zeros(length(a2a),4);
for i = 1:length(a2a)
    firstLastStore(i,1) = DSstore(1,i);
    firstLastStore(i,3) = NSstore(1,i);
    %find the last numeric value
    findNums = find(~isnan(DSstore(:,i)));
    firstLastStore(i,2) = DSstore(findNums(end),i);
    firstLastStore(i,4) = NSstore(findNums(end),i);
end

hFig = figure;
hold on
for i = 1:length(a2a)
    plot(firstLastStore(i,1:2),'b.-')
    plot(firstLastStore(i,3:4),'r.-')
end
errorbar(mean(firstLastStore(:,1:2)),std(firstLastStore(:,1:2))/sqrt(length(a2a)),'b.-','LineWidth',2)
errorbar(mean(firstLastStore(:,3:4)),std(firstLastStore(:,3:4))/sqrt(length(a2a)),'r.-','LineWidth',2)
xlim([0 3])
xlabel('First vs Last Day')
ylabel('Photometry Signal (Zscored)')
title('A2A First and Last Day')
savefig(hFig,'a2afirstVSlastPhotValue');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'a2afirstVSlastPhotValue','-dpdf','-r0')

%now we can do ratio of change between the two
firstLastRatio = zeros(length(a2a),2);
firstLastRatio(:,1) = firstLastStore(:,2)./firstLastStore(:,1);
firstLastRatio(:,2) = firstLastStore(:,4)./firstLastStore(:,3);

hFig = figure;
hold on
for i = 1:length(a2a)
    plot(firstLastRatio(i,:),'k')
end
plot([ones(length(a2a),1)],firstLastRatio(:,1),'b.')
plot([ones(length(a2a),1)*2],firstLastRatio(:,2),'r.')
xlim([0 3])
xlabel('DS (1) and NS(2)')
ylabel('Last/First Day Photometry')
title('A2A Ratio Change First and Last Day')
savefig(hFig,'a2afirstVSlastPhotRatio');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'a2afirstVSlastPhotRatio','-dpdf','-r0')

%we can also try linear regression fit for all the sets of behavioral data
slopeStore = [];
for i = 1:length(a2a)
    x = find(~isnan(DSstore(:,i)));
    y1 = DSstore(x,i);
    P = polyfit(x,y1,1);
    slopeStore(i,1) = P(1);
    y1 = NSstore(x,i);
    P = polyfit(x,y1,1);
    slopeStore(i,2) = P(1);
end


hFig = figure;
hold on
for i = 1:length(a2a)
    plot(slopeStore(i,:),'k')
end
plot([ones(length(a2a),1)],slopeStore(:,1),'b.')
plot([ones(length(a2a),1)*2],slopeStore(:,2),'r.')
xlim([0 3])
xlabel('DS (1) and NS(2)')
ylabel('Slope for Overall Behavior')
title('A2A Slope First and Last Day')
savefig(hFig,'a2afirstVSlastPhotSlope');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'a2afirstVSlastPhotSlope','-dpdf','-r0')


%% okay, now lets go after the auc of photometry predicting trial type!

%this is teh fifth column

hFig = figure
hold on
meanStore = nan(14,length(d1));
for i = 1:length(d1)
    plot(fullStore{d1(i),5})
    lengthFind = length(fullStore{d1(i),5});
    meanStore(1:lengthFind,i) = fullStore{d1(i),5};
end
plot(nanmean(meanStore'),'b','LineWidth',2)
xlabel('Days of Training')
ylabel('AUC Score for Photometry Predicting Trial Type')
title('AUC for D1 Over Training')
savefig(hFig,'d1PhotAUC');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'d1PhotAUC','-dpdf','-r0')


hFig = figure
hold on
meanStore = nan(14,length(d1));
for i = 1:length(d1)
    plot(fullStore{d1(i),5})
    lengthFind = length(fullStore{d1(i),5});
    meanStore(1:lengthFind,i) = fullStore{d1(i),5};
end
plot(nanmean(meanStore'),'b','LineWidth',2)
xlabel('Days of Training')
ylabel('AUC Score for Photometry Predicting Trial Type')
title('AUC for D1 Over Training')
savefig(hFig,'d1PhotAUC');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'d1PhotAUC','-dpdf','-r0')

hFig = figure
hold on
meanStore = nan(14,length(a2a));
for i = 1:length(a2a)
    plot(fullStore{a2a(i),5})
    lengthFind = length(fullStore{a2a(i),5});
    meanStore(1:lengthFind,i) = fullStore{a2a(i),5};
end
plot(nanmean(meanStore'),'b','LineWidth',2)
xlabel('Days of Training')
ylabel('AUC Score for Photometry Predicting Trial Type')
title('AUC for A2A Over Training')
savefig(hFig,'a2aPhotAUC');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'a2aPhotAUC','-dpdf','-r0')


%% eliminate mice that dont learn and plot this as well. 
d1Fix = d1;
d1Fix(5) = [];
a2aFix = a2a;
a2aFix(7) = [];

hFig = figure
hold on
meanStore = nan(14,length(d1Fix));
for i = 1:length(d1Fix)
    plot(fullStore{d1Fix(i),5})
    lengthFind = length(fullStore{d1Fix(i),5});
    meanStore(1:lengthFind,i) = fullStore{d1Fix(i),5};
end
plot(nanmean(meanStore'),'b','LineWidth',2)
xlabel('Days of Training')
ylabel('AUC Score for Photometry Predicting Trial Type')
title('AUC for D1 Over Training')
savefig(hFig,'d1FixPhotAUC');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'d1FixPhotAUC','-dpdf','-r0')

hFig = figure
hold on
meanStore = nan(14,length(a2aFix));
for i = 1:length(a2aFix)
    plot(fullStore{a2aFix(i),5})
    lengthFind = length(fullStore{a2aFix(i),5});
    meanStore(1:lengthFind,i) = fullStore{a2aFix(i),5};
end
plot(nanmean(meanStore'),'b','LineWidth',2)
xlabel('Days of Training')
ylabel('AUC Score for Photometry Predicting Trial Type')
title('AUC for A2A Over Training')
savefig(hFig,'a2aFixPhotAUC');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'a2aFixPhotAUC','-dpdf','-r0')



