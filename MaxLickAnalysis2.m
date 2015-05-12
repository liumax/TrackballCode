
clear
clc

[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select Lick Data')
addpath(PathName);
load(FileName)

upperbound=4; %acceptable stds above mean for lick duration 
lowerbound=2; %acceptable stds below mean for lick duration
waterthresh=1; %threshold for TTL signal from solenoid
cuethresh=1;
window=[-5,5]; %window around sound onset for raster

freqbin=0.2; %bin size for calculating frequency of licking, used for psth

%MAY WANT TO TWEAK LICK THRESHOLD. 2 IS DEFINITELY ON LOW SIDE
LickThreshold=2.0;       %Voltage Threshold for Lickometer

%Imports lick data
LickVoltages=LickData(:,1);
smVolt=smooth(LickVoltages,10);


%This is to check all values to make sure randomization is done
%appropriately

figure
subplot(2,1,1)
hist(master(:,3),20)
title('ITI times')
xlabel('Time (s)')
ylabel('Trials')
subplot(2,1,2)
hist(master(:,5)-master(:,4),20)
title('Time between Cue and Water')
xlabel('Time (s)')
ylabel('Trials')

%Convert smoothLick to binary. Visually confirmed that 2 is a good cutoff
binLick=zeros(length(smVolt),1);
binLick(smVolt>LickThreshold)=1;
binLick(smVolt<LickThreshold)=0;

%Finds lick onsets and offsets
startind=find(diff(binLick)>0);
endind=find(diff(binLick)<0);

%Checks to make sure that all licks are completed
indcheck=length(startind)-length(endind)

%Throws error if licks are not completed
if abs(indcheck)>1
    disp 'Error in indices, check startind and endind'
    pause
end

%Automatically deletes one index to correct error.
if indcheck<0
    endind(end)=[];
elseif indcheck>0
    startind(end)=[];
end

%calculation of lick duration and fitting of gaussian
lickdur=endind-startind;
[y,x]=hist(lickdur,1000);
q=fit(x',y','gauss1');
lickmean=q.b1;
lickstd=q.c1;

%setting thresholds of licking duration, based on calculated values AND
%preset stdev limits
lickhi=lickmean+upperbound*lickstd;
licklo=lickmean-lowerbound*lickstd;

%finds all bad licks, then removes them
badind=lickdur>lickhi | lickdur<licklo;
badlicks=find(badind==1);
for q=1:length(badlicks)
    binLick(startind(badlicks(q)):endind(badlicks(q)),1)=0;
end

%converts lickstart into double. Lickstart is used as time of licking
lickstart=diff(binLick)>0;
lickstart(end+1)=0;
lickstart=+lickstart;

%imports time data
licktimes=LickTime(startind);

%calculates interlick period
interlick=diff(licktimes);
lickfreq=1./interlick;

%calculates median and frequencies
lickmedian=median(lickfreq);
lickmax=max(lickfreq);
[y,x]=hist(lickfreq,1000);

%plots figures
figure
subplot(2,1,1)
hist(lickdur,1000)
title('Lick Durations')
subplot(2,1,2)
hist(lickfreq,1000);
title(['Inter-Lick Interval in Hz with Median of ' ,num2str(lickmedian), 'Hz'])


%%Finds Start of Water Delivery
binWater=LickData(:,3);
binWater(binWater>waterthresh)=1;
binWater(binWater<waterthresh)=0;
startWater=find(diff(binWater)==1);

%%Finds Start of Cue
binCue=LickData(:,2);
binCue(binCue>cuethresh)=1;
binCue(binCue<cuethresh)=0;
startCue=find(diff(binCue)==1);

rasteraxis=[window(1):0.001:window(2)-0.001];

%Makes Raster of all Rewarded Trials Relative to Reward Delivery
waterprep=zeros(1000*(abs(window(1))+abs(window(2))),length(startWater));
%makes x by y array of licking data around water starts, with x being
%number of samples (determined by window) and y being number of trials.
for i=1:length(startWater)
    waterprep(:,i)=lickstart(startWater(i)+(window(1)*1000):startWater(i)+(window(2)*1000)-1);
end
%generates raster
counter=1;
for i=1:size(waterprep,2)
    findlicks=[];
    findlicks=find(waterprep(:,i)==1);
    findlicksize=size(findlicks,1);
    waterraster(counter:counter+findlicksize-1,2)=findlicks;
    waterraster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end
%replaces raster index with times, so can graph, also calculates histogram
%to bin licks. 
waterraster(:,2)=rasteraxis(waterraster(:,2));
[watercounts,watercenters]=hist(waterraster(:,2),(window(2)-window(1))/freqbin);


%Makes Raster of all Trials Relative to Cue Delivery
cueprep=zeros(1000*(abs(window(1))+abs(window(2))),length(startCue));

%makes x by y array of licking data around cue starts, with x being
%number of samples (determined by window) and y being number of trials.
for i=1:length(startCue)
    cueprep(:,i)=lickstart(startCue(i)+(window(1)*1000):startCue(i)+(window(2)*1000)-1);
end
%generates raster
counter=1;
for i=1:size(cueprep,2)
    findlicks=[];
    findlicks=find(cueprep(:,i)==1);
    findlicksize=size(findlicks,1);
    cueraster(counter:counter+findlicksize-1,2)=findlicks;
    cueraster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end

%replaces raster index with times, so can graph, also calculates histogram
%to bin licks. 
cueraster(:,2)=rasteraxis(cueraster(:,2));
[cuecounts,cuecenters]=hist(cueraster(:,2),(window(2)-window(1))/freqbin);

figure
subplot(2,2,1)
plot(waterraster(:,2),waterraster(:,1),'b.')
title('Raster Relative to Water')
xlabel('Time (s)')
ylabel('Rasters')

subplot(2,2,2)
plot(cueraster(:,2),cueraster(:,1),'r.')
title('Raster Relative to Cue')
xlabel('Time (s)')
ylabel('Rasters')

subplot(2,2,3)
plot(watercenters,watercounts,'b')
title('Total Licks Relative to Water')
xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
ylabel('Number of Licks')

subplot(2,2,4)
plot(cuecenters,cuecounts,'r')
title('Total Licks Relative to Cue')
xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
ylabel('Number of Licks')

%%%Rasters of Big vs Small Trials
rewardwaters=master(:,1);
rewardwaters(master(:,6)==0)=[];

bigrewards=find(rewardwaters==max(rewardwaters));
smallrewards=find(rewardwaters==min(rewardwaters));

%Makes Raster for Big/Small in Two Colors Relative to Water

bwprep=zeros(size(waterprep,1),size(bigrewards,1));
swprep=zeros(size(waterprep,1),size(smallrewards,1));

for i=1:length(bigrewards)
    bwprep(:,i)=waterprep(:,(bigrewards(i)));
end

counter=1;
for i=1:size(bwprep,2)
    findlicks=[];
    findlicks=find(bwprep(:,i)==1);
    findlicksize=size(findlicks,1);
    bwraster(counter:counter+findlicksize-1,2)=findlicks;
    bwraster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end
bwraster(:,2)=rasteraxis(bwraster(:,2));
bwraster(:,1)=bwraster(:,1)+size(master,1)/2;

[bwcounts,bwcenters]=hist(bwraster(:,2),(window(2)-window(1))/freqbin);

for i=1:length(smallrewards)
    swprep(:,i)=waterprep(:,(smallrewards(i)));
end

counter=1;
for i=1:size(swprep,2)
    findlicks=[];
    findlicks=find(swprep(:,i)==1);
    findlicksize=size(findlicks,1);
    swraster(counter:counter+findlicksize-1,2)=findlicks;
    swraster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end
swraster(:,2)=rasteraxis(swraster(:,2));

[swcounts,swcenters]=hist(swraster(:,2),(window(2)-window(1))/freqbin);

%Makes Raster for Big/Small in Two Colors Relative to Cue
bigcue=find(master(:,6)==1&master(:,1)==max(master(:,1)));
smallcue=find(master(:,6)==1&master(:,1)==min(master(:,1)));

bcprep=zeros(size(cueprep,1),size(bigcue,1));
scprep=zeros(size(cueprep,1),size(smallcue,1));

for i=1:length(bigcue)
    bcprep(:,i)=cueprep(:,(bigcue(i)));
end

counter=1;
for i=1:size(bcprep,2)
    findlicks=[];
    findlicks=find(bcprep(:,i)==1);
    findlicksize=size(findlicks,1);
    bcraster(counter:counter+findlicksize-1,2)=findlicks;
    bcraster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end
bcraster(:,2)=rasteraxis(bcraster(:,2));
bcraster(:,1)=bcraster(:,1)+size(master,1)/2;

[bccounts,bccenters]=hist(bcraster(:,2),(window(2)-window(1))/freqbin);

for i=1:length(smallcue)
    scprep(:,i)=cueprep(:,(smallcue(i)));
end

counter=1;
for i=1:size(scprep,2)
    findlicks=[];
    findlicks=find(scprep(:,i)==1);
    findlicksize=size(findlicks,1);
    scraster(counter:counter+findlicksize-1,2)=findlicks;
    scraster(counter:counter+findlicksize-1,1)=i;
    counter=counter+findlicksize;
end
scraster(:,2)=rasteraxis(scraster(:,2));

[sccounts,sccenters]=hist(scraster(:,2),(window(2)-window(1))/freqbin);

%Plots data out.
figure
subplot(2,2,1)
plot(swraster(:,2),swraster(:,1),'k.',bwraster(:,2),bwraster(:,1),'m.')
title('Rasters of Small (black) and Large (Magenta) Rewards Relative to Water')
xlabel('Time (s)')
ylabel('Rasters')
subplot(2,2,3)
plot(swcenters,swcounts,'k',bwcenters,bwcounts,'m')
title('Histograms of Licks Relative to Water ')
xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
ylabel('Number of Licks')
subplot(2,2,2)
plot(scraster(:,2),scraster(:,1),'k.',bcraster(:,2),bcraster(:,1),'m.')
title('Rasters of Small (black) and Large (Magenta) Rewards Relative to Cue')
xlabel('Time (s)')
ylabel('Rasters')
subplot(2,2,4)
plot(sccenters,sccounts,'k',bccenters,bccounts,'m')
title('Histograms of Licks Relative to Cue ')
xlabel(['Time (s) in ',num2str(freqbin),'sec bins'])
ylabel('Number of Licks')


%%Calculates total licks between cue and water delivery (excludes
%%consummatory licks)

antilicks=zeros(size(startCue,1),2);
antilicksb=zeros(size(bigcue,1),1);
antilickss=zeros(size(smallcue,1),1);

for i=1:size(startCue,1)
    antilicks(i,1)=sum(lickstart(startCue(i):startWater(i)));
    antilicks(i,2)=master(i,1);
end

for i=1:size(antilicksb,1)
    antilicksb(i)=sum(lickstart(startCue(bigcue(i)):startWater(bigrewards(i))));
end

for i=1:size(antilickss,1)
    antilickss(i)=sum(lickstart(startCue(smallcue(i)):startWater(smallrewards(i))));
end

figure
subplot(2,1,1)
hist(antilickss,10)
title(['Histogram of Anticipatory Licks, Small Reward, Mean of ',num2str(mean(antilickss)),' Licks/Trial'])
xlabel('Number of Licks')
subplot(2,1,2)
hist(antilicksb,10)
title(['Histogram of Anticipatory Licks, Big Reward, Mean of ',num2str(mean(antilicksb)),' Licks/Trial'])
xlabel('Number of Licks')

%find latency to first lick
licklatb=cell(size(bigcue,1),1);
licklats=cell(size(bigcue,1),1);

for i=1:size(bigcue,1)
    licklatb{i}=find(lickstart(startCue(bigcue(i)):startWater(bigrewards(i)))==1,1,'first');
end

for i=1:size(smallcue,1)
    licklats{i}=find(lickstart(startCue(smallcue(i)):startWater(smallrewards(i)))==1,1,'first');
end

licklatb=cell2mat(licklatb(~cellfun('isempty',licklatb)));
licklats=cell2mat(licklats(~cellfun('isempty',licklats)));
licklatbmean=mean(licklatb);
licklatsmean=mean(licklats);

figure
subplot(2,1,1)
hist(licklats)
title(['Histogram of Lick Latency, Small Rewards Mean of ' num2str(licklatsmean) ' msec'])
xlabel('Milliseconds')
subplot(2,1,2)
hist(licklatb)
title(['Histogram of Lick Latency, Big Rewards Mean of ' num2str(licklatbmean) ' msec'])
xlabel('Milliseconds')






