%THis is code to make a fake spike train!!

load('181205_ML180928G_R_AudStr_pen2_rec2_3mWPVHaloTuningAltLaserFullTuningAnalysis.mat')

%these are our targets. 
figure
plot(s.nt3cluster1.AverageWaveForms)
%looking at channel 1. 
figure
plot(s.nt5cluster1.AverageWaveForms)
%looking at channel 3.


%interpolate really fine. 1:10

for i = 1:4
    interpFSI(:,i) = interp1([1:1:40],s.nt3cluster1.AverageWaveForms(:,i),[1:0.1:40]);
    interpMSN(:,i) = interp1([1:1:40],s.nt5cluster1.AverageWaveForms(:,i),[1:0.1:40]);
end

interpFSI = -1*interpFSI;
interpMSN = -1*interpMSN;

%generate fake train. waves are now 390 samples long, which corresponds to
%1.3 ms. 


tester = zeros(391*20,4);
fsiTimes = [400 1600 2100 3000 3400 4800 6800 7300];
msnTimes = [800 3200 5000];

for i = 1:length(fsiTimes)
    tester(fsiTimes(i):fsiTimes(i) +391 -1,:) = tester(fsiTimes(i):fsiTimes(i) +391 -1,:) + interpFSI;
end

for i = 1:length(msnTimes)
    tester(msnTimes(i):msnTimes(i) +391 -1,:) = tester(msnTimes(i):msnTimes(i) +391 -1,:) + interpMSN;
end

hFig = figure;
hold on
for i = 1:4
    plot(tester(:,i) - 100*i,'k','LineWidth',2)
end

set(gca,'TickDir','out')
spikeGraphName = 'ExampleSpikeTrain';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
hold on
for i = 1:4
    plot(tester(:,i) - 100*i,'Color',[0.7 0.7 0.7],'LineWidth',2)
    if i == 1
        for j = 1:length(fsiTimes)
            plot([round(fsiTimes(j)):round(fsiTimes(j) + 391 -1)],tester([round(fsiTimes(j)):round(fsiTimes(j) + 391 -1)],i)- 100*i,'r','LineWidth',2)
        end
    elseif i > 1 & i < 4
        for j = 1:length(fsiTimes)
            plot([round(fsiTimes(j)):round(fsiTimes(j) + 391 -1)],tester([round(fsiTimes(j)):round(fsiTimes(j) + 391 -1)],i)- 100*i,'r','LineWidth',2)
        end
        for j = 1:length(msnTimes)
            plot([round(msnTimes(j)):round(msnTimes(j) + 391 -1)],tester([round(msnTimes(j)):round(msnTimes(j) + 391 -1)],i)- 100*i,'k','LineWidth',2)
        end
    elseif i == 4
        for j = 1:length(msnTimes)
            plot([round(msnTimes(j)):round(msnTimes(j) + 391 -1)],tester([round(msnTimes(j)):round(msnTimes(j) + 391 -1)],i)- 100*i,'k','LineWidth',2)
        end
    end
end

set(gca,'TickDir','out')
spikeGraphName = 'ExampleSpikeTrainLabeled';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
subplot(1,2,1)
plot(interpFSI(:,1),'r','LineWidth',2)
ylim([-150 100])
set(gca,'TickDir','out')
subplot(1,2,2)
plot(interpMSN(:,3),'k','LineWidth',2)
ylim([-150 100])
set(gca,'TickDir','out')

spikeGraphName = 'ExampleSpikeIndivWave';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


hFig = figure;
anchorFSI = [];
anchorMSN = [];


tester = zeros(300,2);
tester(:,1) = 90;
tester(:,2) = 30;
X = randn([300,1])*12; %seems to max out around 3. 
tester(:,1) = tester(:,1) + X;
X = randn([300,1])*12; %seems to max out around 3. 
tester(:,2) = tester(:,2) + X;
cloudFSI = tester;

tester = zeros(300,2);
tester(:,1) = 30;
tester(:,2) = 120;
X = randn([300,1])*18; %seems to max out around 3. 
tester(:,1) = tester(:,1) + X;
X = randn([300,1])*18; %seems to max out around 3. 
tester(:,2) = tester(:,2) + X;
cloudMSN = tester;

hFig = figure;
hold on
plot(cloudFSI(:,1),cloudFSI(:,2),'k.')
plot(cloudMSN(:,1),cloudMSN(:,2),'k.')
set(gca,'TickDir','out')

spikeGraphName = 'Spike Clusters';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
hold on
plot(cloudFSI(:,1),cloudFSI(:,2),'r.')
plot(cloudMSN(:,1),cloudMSN(:,2),'k.')
set(gca,'TickDir','out')

spikeGraphName = 'Spike Clusters Colored';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%% Now make a bunch of sine waves. 

f=4;
Amp=1;
ts=1/8000;
T=0.25;
t=0:ts:T;
y=sin(2*pi*f*t);
% figure;
% plot(t,y)

sineStore = zeros(length(t),5,16);

for i = 1:5
    amp = i/10;
    for j = 1:16
        f = 4 + 1*(j-1);
        sineStore(:,i,j) = amp*sin(2*pi*f*t);
    end
end

plotStore = zeros(length(t)*(16+4),5);
counter = 1;
for i = 1:16
    plotStore(counter:counter + length(t) - 1,:) = sineStore(:,:,i);
    counter = counter + length(t);
    plotStore(counter:counter + floor(length(t)/4),:) = NaN;
    counter = counter + floor(length(t)/4);
end


hFig = figure;
hold on
for i = 1:5
    plot(plotStore(:,i) + 2*i)
end


spikeGraphName = 'SineWaveSounds';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

