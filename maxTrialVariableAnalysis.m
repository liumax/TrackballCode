function [cleanLicks] = maxTrialVariableAnalysis(trialStates,trialParams,portStates);

lickUpper=150;
lickLower=5;

window=[-5,5]; %window for raster

freqbin=0.2; %bin size for calculating frequency of licking, used for psth

%This is to check all values to make sure randomization is done
%appropriately

% figure
% subplot(2,1,1)
% hist(diff(trialStates.soundOn)/1000,100)
% title('ITI times')
% axis([0 20 0 400])
% xlabel('Time (s)')
% ylabel('Trials')
% subplot(2,1,2)
% hist(trialStates.rewDelivery - trialStates.soundOn,100)
% title('Time between Cue and Water')
% xlabel('Time (s)')
% ylabel('Trials')

lickArray = [portStates.tStamps',portStates.inStates(:,3)];

dupCounter = 1; %this is a placeholder to have a counter update

%This finds all duplicate values in adjacent rows.
for i = 2: size(lickArray,1)
    checker = lickArray(i-1,2) == lickArray(i,2);
    if checker == 1
        dupes(dupCounter) = i;
        dupCounter = dupCounter + 1;
    end
end
%This deletes all duplicate values.
for i = size(dupes,2):-1:1
    lickArray(dupes(i),:) = [];
end
%since all duplicates are gone, values of 1 indicate lick starts, 0
%indicates lick end
lickStart=find(lickArray(:,2) == 1);
lickEnd=find(lickArray(:,2) == 0);
lickDur = lickArray(lickEnd,1) - lickArray(lickStart,1);
badlicks = find(lickDur > lickUpper | lickDur < lickLower);

cleanLicks=lickArray;

for i = size(badlicks,1):-1:1
    cleanLicks(lickStart(badlicks(i)):lickEnd(badlicks(i)),:)=[];
end

%cleanlicks is now licks with no licks that do not meet criteria

cleanStart=find(cleanLicks(:,2) == 1);
cleanEnd=find(cleanLicks(:,2) == 0);

lickPeriod=diff(cleanLicks(cleanStart));
lickFreq=1./(lickPeriod/1000);


%plots figures
% figure
% subplot(2,1,1)
% hist(lickDur,1000)
% title('Lick Durations')
% subplot(2,1,2)
% hist(lickFreq,1000);
% title(['Inter-Lick Interval in Hz'])


%%Finds Start of Water Delivery
startWater=trialStates.rewDelivery;
startWater(startWater<1) = []; %This gets rid of hanger-ons for last trial displayed but not executed

%%Finds Start of Cue
startCue=trialStates.soundOn;
startCue(startCue<1) = [];%This gets rid of hanger-ons for last trial displayed but not executed

startLick = cleanLicks(cleanStart,1);

rasteraxis=[window(1):0.001:window(2)-0.001];

waterPrepper=cell(size(trialStates.rewLength,2),1);%generates open cell array

for i = 1:size(waterPrepper,1)
    waterPrepper{i}=startLick - startWater(i);%subtracts water point from all times for licks
    placer=find(abs(waterPrepper{i})<5000);%finds licks within window of 5sec
    waterPrepper{i}=waterPrepper{i}(placer);%eliminates other licks
end

counter = 1;
for i=1:size(waterPrepper,1)
    findlicksize=size(waterPrepper{i},1);
    waterRaster(counter:counter+findlicksize-1,2)=waterPrepper{i};
    waterRaster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end

[waterCounts,waterCenters]=hist(waterRaster(:,2),(window(2)-window(1))/freqbin);


cuePrepper=cell(size(trialStates.rewLength,2),1);

for i = 1:size(cuePrepper,1)
    cuePrepper{i}=startLick - startCue(i);
    placer=find(abs(cuePrepper{i})<5000);
    cuePrepper{i}=cuePrepper{i}(placer);
end

counter = 1;
for i=1:size(cuePrepper,1)
    findlicksize=size(cuePrepper{i},1);
    cueRaster(counter:counter+findlicksize-1,2)=cuePrepper{i};
    cueRaster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end

[cueCounts,cueCenters]=hist(cueRaster(:,2),(window(2)-window(1))/freqbin);

% % figure
% % subplot(2,2,1)
% % plot(waterRaster(:,2),waterRaster(:,1),'b.')
% % title('Raster Relative to Water')
% % xlabel('Time (s)')
% % ylabel('Rasters')
% % 
% % subplot(2,2,2)
% % plot(cueRaster(:,2),cueRaster(:,1),'r.')
% % title('Raster Relative to Cue')
% % xlabel('Time (s)')
% % ylabel('Rasters')
% % 
% % subplot(2,2,3)
% % plot(waterCenters,waterCounts,'b')
% % title('Total Licks Relative to Water')
% % xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
% % ylabel('Number of Licks')
% % 
% % subplot(2,2,4)
% % plot(cueCenters,cueCounts,'r')
% % title('Total Licks Relative to Cue')
% % xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
% % ylabel('Number of Licks')

rewardSize=trialStates.rewLength;
minSize=str2double(trialParams.minRew);
maxSize=str2double(trialParams.maxRew);

bigRewards=find(rewardSize == maxSize);
smallRewards=find(rewardSize == minSize);

%below calculates rasters relative to water for big and small
bwprep=cell(size(bigRewards,1),1);
swprep=cell(size(smallRewards,1),1);

for i=1:length(bigRewards)
    bwprep{i}=waterPrepper{bigRewards(i)};
end

counter=1;
for i=1:size(bwprep,2)
    findlicksize=size(bwprep{i},1);
    bwRaster(counter:counter+findlicksize-1,2)=bwprep{i};
    bwRaster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end

[bwCounts,bwCenters]=hist(bwRaster(:,2),(window(2)-window(1))/freqbin);


for i=1:length(smallRewards)
    swprep{i}=waterPrepper{smallRewards(i)};
end

counter=1;
for i=1:size(swprep,2)
    findlicksize=size(swprep{i},1);
    swRaster(counter:counter+findlicksize-1,2)=swprep{i};
    swRaster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end

[swCounts,swCenters]=hist(swRaster(:,2),(window(2)-window(1))/freqbin);

%below calculates rasters and histograms relative to cue
bcprep=cell(size(bigRewards,1),1);
scprep=cell(size(smallRewards,1),1);

for i=1:length(bigRewards)
    bcprep{i}=cuePrepper{bigRewards(i)};
end

counter=1;
for i=1:size(bcprep,2)
    findlicksize=size(bcprep{i},1);
    bcRaster(counter:counter+findlicksize-1,2)=bcprep{i};
    bcRaster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end

[bcCounts,bcCenters]=hist(bcRaster(:,2),(window(2)-window(1))/freqbin);


for i=1:length(smallRewards)
    scprep{i}=cuePrepper{smallRewards(i)};
end

counter=1;
for i=1:size(scprep,2)
    findlicksize=size(scprep{i},1);
    scRaster(counter:counter+findlicksize-1,2)=scprep{i};
    scRaster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end

[scCounts,scCenters]=hist(scRaster(:,2),(window(2)-window(1))/freqbin);


%Plots data out.
figure
subplot(2,2,1)
plot(swRaster(:,2),swRaster(:,1),'k.',bwRaster(:,2),bwRaster(:,1),'r.')
title('Rasters of Small (black) and Large (Red) Rewards Relative to Water')
xlabel('Time (s)')
ylabel('Rasters')
subplot(2,2,3)
plot(swCenters,swCounts,'k',bwCenters,bwCounts,'r')
title('Histograms of Licks Relative to Water ')
xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
ylabel('Number of Licks')
subplot(2,2,2)
plot(scRaster(:,2),scRaster(:,1),'k.',bcRaster(:,2),bcRaster(:,1),'r.')
title('Rasters of Small (black) and Large (Red) Rewards Relative to Cue')
xlabel('Time (s)')
ylabel('Rasters')
subplot(2,2,4)
plot(scCenters,scCounts,'k',bcCenters,bcCounts,'r')
title('Histograms of Licks Relative to Cue ')
xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
ylabel('Number of Licks')


figure
hold on
plot(trialStates.preLick,'k','linewidth',2)
plot(trialStates.postLick,'b','linewidth',2)
plot(trialStates.antiLick,'g','linewidth',2)
axis([0 400 0 30])
title('Licks Over Session')
xlabel('Trials')
ylabel('Licks')
hold off

end





