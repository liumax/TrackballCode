
%% Okay, now that we are done with that shit, lets look at group data and produce some basic group figures
%  fullMaster = zeros(850,10,4);
%  fullMaster(:,:,1) = masterI;
%  fullMaster(:,:,2) = masterA;
%  fullMaster(:,:,3) = masterC;
%  fullMaster(:,:,4) = masterD;
%  
%  fullStruct.i = structI;
%  fullStruct.a = structA;
%  fullStruct.c = structC;
%  fullStruct.d = structD;


%so now lets try and determine how long it took animals to learn this
%shit. based on previous graphing, I think the difference in anticipatory
%licks isnt a bad measure. We can go by the day. 

load('CombinedDataset.mat')

for i = 1:4
 dailyTrace(i,:) = (fullSummary(:,2,i)-fullSummary(:,5,i))./(fullSummary(:,2,i)+fullSummary(:,5,i));
end

hFig = figure
plot(dailyTrace(:,2:end)','LineWidth',2)
xlabel('Days of Training')
ylabel('Lick Preference')
title('Learning of Pavlovian Task Over Time')
legend('D1.1','D1.2','A2A.1','A2A.2','Location','northwest')

savefig(hFig,'LickingPreference');
    
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'LickingPreference','-dpdf','-r0')

%average across animals, include standard deviationp

meanDaily = mean(dailyTrace(:,2:end));
steDaily = std(dailyTrace(:,2:end))/sqrt(4);

figure
plot(meanDaily,'LineWidth',3)
hold on
plot(meanDaily + steDaily,'LineWidth',2)
plot(meanDaily - steDaily,'LineWidth',2)
% plot(dailyTrace(:,2:end)')
xlabel('Days of Training')
ylabel('Lick Preference')
title('Lick Preference Over Time')
ylim([-1 1])

%baseline licking?

for i = 1:4
 baselineTrace(i,:) = fullSummary(:,1,i);
end

hFig = figure
plot(baselineTrace(:,2:end)','LineWidth',2)
xlabel('Days of Training')
ylabel('Baseline Licks')
title('Baseline Licking Over Time')
legend('D1.1','D1.2','A2A.1','A2A.2','Location','northwest')

savefig(hFig,'BaselineLicks');
    
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'BaselineLicks','-dpdf','-r0')

meanBase = mean(baselineTrace(:,2:end));
steBase = std(baselineTrace(:,2:end))/2;

figure
plot(meanBase,'LineWidth',3)
hold on
plot(meanBase + steBase,'LineWidth',2)
plot(meanBase - steBase,'LineWidth',2)
plot(baselineTrace(:,2:end)')
xlabel('Days of Training')
ylabel('Baseline Licks')
title('Baseline Licking Over Time')

%licking to wrong cue?

for i = 1:4
 wrongTrace(i,:) = fullSummary(:,5,i);
end

figure
plot(wrongTrace(:,2:end)','LineWidth',2)
xlabel('Days of Training')
ylabel('Difference in Anticipatory Licks')
title('Learning of Pavlovian Task Over Time')


meanWrong = mean(wrongTrace(:,2:end));
steWrong = std(wrongTrace(:,2:end))/2;


figure
plot(meanWrong,'LineWidth',2)
hold on
plot(meanWrong + steWrong)
plot(meanWrong - steWrong)
xlabel('Days of Training')
ylabel('Anti Licks Wrong')
% title('Baseline Licking Over Time')


subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,2,1)
plot(squeeze(sumStore(2:end,2,1)))
hold on
plot(squeeze(sumStore(2:end,5,1)),'r')
ylabel('Anticipatory Licks')
xlabel('Training Days')
legend('DS','NS','Location','northwest')
subplot(2,2,3)
plot(squeeze(sumStore(2:end,2,2)))
hold on
plot(squeeze(sumStore(2:end,5,2)),'r')
ylabel('Anticipatory Licks')
xlabel('Training Days')
legend('DS','NS','Location','northwest')
subplot(2,2,2)
plot(squeeze(sumStore(2:end,2,3)))
hold on
plot(squeeze(sumStore(2:end,5,3)),'r')
ylabel('Anticipatory Licks')
xlabel('Training Days')
legend('DS','NS','Location','northwest')
subplot(2,2,4)
plot(squeeze(sumStore(2:end,2,4)))
hold on
plot(squeeze(sumStore(2:end,5,4)),'r')
ylabel('Anticipatory Licks')
xlabel('Training Days')
legend('DS','NS','Location','northwest')
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'anticipatoryLicksOverall','-djpeg90','-r0')

































