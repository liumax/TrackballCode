
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

figure
plot(dailyTrace(:,2:end)','LineWidth',2)
xlabel('Days of Training')
ylabel('Difference in Anticipatory Licks')
title('Learning of Pavlovian Task Over Time')

%average across animals, include standard deviationp

meanDaily = mean(dailyTrace(:,2:end));
steDaily = std(dailyTrace(:,2:end))/sqrt(4);

figure
plot(meanDaily,'LineWidth',3)
hold on
plot(meanDaily + steDaily,'LineWidth',2)
plot(meanDaily - steDaily,'LineWidth',2)
plot(dailyTrace(:,2:end)')
xlabel('Days of Training')
ylabel('Lick Preference')
title('Lick Preference Over Time')
ylim([-1 1])

%baseline licking?

for i = 1:4
 baselineTrace(i,:) = fullSummary(:,1,i);
end

figure
plot(baselineTrace(:,2:end)','LineWidth',2)
xlabel('Days of Training')
ylabel('Baseline Licks')
title('Baseline Licking Over Time')

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

for i = 1:4
 trialTrace(i,:) = smooth(fullMaster(:,2,i)-fullMaster(:,5,i),25);
end

figure
plot(trialTrace')
hold on
plot([50 50],[-2 5],'r')
plot([150 150],[-2 5],'r')
plot([250 250],[-2 5],'r')
plot([350 350],[-2 5],'r')
plot([450 450],[-2 5],'r')
plot([550 550],[-2 5],'r')
plot([650 650],[-2 5],'r')
plot([750 750],[-2 5],'r')
plot([850 850],[-2 5],'r')


%now lets look at basic data from taskNoRew: what do we see?



tester = fullStruct.i.HiRaster(:,1:50);



































