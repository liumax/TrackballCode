
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
     dailyTrace(i,:) = fullSummary(:,2,i)-fullSummary(:,5,i);
 end
 
 figure
 plot(dailyTrace(:,2:end)')
 xlabel('Days of Training')
 ylabel('Difference in Anticipatory Licks')
 title('Learning of Pavlovian Task Over Time')
 
 %baseline licking?
 
 for i = 1:4
     baselineTrace(i,:) = fullSummary(:,1,i);
 end
 
 figure
 plot(baselineTrace')
 xlabel('Days of Training')
 ylabel('Baseline Licks')
 title('Baseline Licking Over Time')
 
 
 for i = 1:4
     trialTrace(i,:) = smooth(fullMaster(:,2,i)-fullMaster(:,5,i),31);
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

 
 