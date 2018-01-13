%this is test code to pull out 50 pulse laser high low data

fileNames = what;
fileNames = fileNames.mat;
counter = 1;
for i = 1:length(fileNames)
    load(fileNames{i})
    
    %pull master sheet
    masterDesig = s.MasterDesigs;
    masterSheet(counter:counter+size(s.MasterSheet,1)-1,:) = s.MasterSheet;
    
    %pull average rates
    for j = 1:length(s.DesignationName)
        averageHistLow(counter-1+j,:) = s.(s.DesignationName{j}).LaserHistLow;
        averageHistHi(counter-1+j,:) = s.(s.DesignationName{j}).LaserHistHi;
        averageZLow(counter-1+j,:) = (s.(s.DesignationName{j}).LaserHistLow - mean(s.(s.DesignationName{j}).HistogramLaser(1:100)))/mean(s.(s.DesignationName{j}).HistogramLaser(1:100));
        averageZHi(counter-1+j,:) = (s.(s.DesignationName{j}).LaserHistHi - mean(s.(s.DesignationName{j}).HistogramLaser(1:100)))/mean(s.(s.DesignationName{j}).HistogramLaser(1:100));
    end
    counter = counter + size(s.MasterSheet,1);
end


pvs = find(masterSheet(:,1) == 1);
msns = find(masterSheet(:,1) == 0);

msnLowRateMean = mean(averageHistLow(msns,:));
msnHiRateMean = mean(averageHistHi(msns,:));

msnLowZMean = mean(averageZLow(msns,:));
msnHiZMean = mean(averageZHi(msns,:));

figure
hold on
plot(msnLowRateMean,'b','LineWidth',2)
plot(msnHiRateMean,'g','LineWidth',2)
figure
hold on
plot(msnLowZMean,'b','LineWidth',2)
plot(msnHiZMean,'g','LineWidth',2)


%calculate modulation index
modLow = (masterSheet(:,20) - (mean(masterSheet(:,18:19)'))')./(masterSheet(:,20) + (mean(masterSheet(:,18:19)'))');
modHi = (masterSheet(:,23) - (mean(masterSheet(:,21:22)'))')./(masterSheet(:,23) + (mean(masterSheet(:,21:22)'))');

figure
subplot(2,1,1)
hist(modLow(msns),[-1:0.1:1])
xlim([-1 1])
title('Modulation Index MSNs, Low Power (3 mW)')
subplot(2,1,2)
hist(modHi(msns),[-1:0.1:1])
xlim([-1 1])
title('Modulation Index MSNs, Hi Power (15 mW)')

figure
subplot(2,1,1)
hist(modLow(pvs),[-1:0.1:1])
xlim([-1 1])
title('Modulation Index PV, Low Power (3 mW)')
subplot(2,1,2)
hist(modHi(pvs),[-1:0.1:1])
xlim([-1 1])
title('Modulation Index PV, Hi Power (15 mW)')
