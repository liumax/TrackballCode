cd Z:\Max\180516EPAnalysis\loco\180612AllAnalysisFiles

load('180514_ML180501A_L_EP_4825_multiStimSecondLightToneAnalysis.mat')
bigSaver = masterData;
bigSaver(:,13) = 1;
load('180522_ML180501D_L_EP_4861_justLocoLightToneAnalysis.mat')
bigSaver(38:64,1:12) = masterData;
bigSaver(38:64,13) = 2;
load('180523_ML180501C_L_EP_4890_firstLocoLightToneAnalysis.mat')
bigSaver(65:94,1:12) = masterData;
bigSaver(65:94,13) = 3;
load('180523_ML180501C_R_EP_5020_firstLocoLightToneAnalysis.mat')
bigSaver(95:114,1:12) = masterData;
bigSaver(95:114,13) = 4;
load('180524_ML180501E_R_EP_4942_locoFirstLightToneAnalysis.mat')
bigSaver(115:128,1:12) = masterData;
bigSaver(115:128,13) = 5;

epFind = find(bigSaver(:,12) == 1);
otherFind = find(bigSaver(:,12) == 0);

figure
plot(bigSaver(epFind,2),bigSaver(epFind,3),'r.')
hold on
plot(bigSaver(otherFind,2),bigSaver(otherFind,3),'k.')

%look at locomotion AUC
figure
subplot(2,1,1)
hist(bigSaver(epFind,10),[0:0.05:1])
xlim([0 1])
title('Putative EP')
subplot(2,1,2)
hist(bigSaver(otherFind,10),[0:0.05:1])
xlim([0 1])
title('Putative Thal')

%look at modulation at locomotor start
figure
subplot(2,1,1)
hist(bigSaver(epFind,7),[-1:0.05:1])
xlim([-1 1])
title('Putative EP')
subplot(2,1,2)
hist(bigSaver(otherFind,7),[-1:0.05:1])
xlim([-1 1])
title('Putative Thal')

%look at loco ends. 

figure
subplot(2,1,1)
hist(bigSaver(epFind,9),[-1:0.05:1])
xlim([-1 1])
title('Putative EP')
subplot(2,1,2)
hist(bigSaver(otherFind,9),[-1:0.05:1])
xlim([-1 1])
title('Putative Thal')







