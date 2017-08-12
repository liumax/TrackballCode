

%This code is meant to extract further information from two tone pavlov
%basic analysis

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);



testNames = what;
testNames = testNames.mat;

bigMaster = zeros(100,10);
sumMaster = zeros(5,2);
bigStruct = struct;
masterInd = 1;

trialNumStore = zeros(5,2);


for bigInd = 1:length(testNames)
    load(testNames{bigInd})
    %first, lets get all of our time points.

    %get ITI, tone times, reward delivery times

    timeITI = s.MBED.ToneDelivery;
    timeHi = s.MBED.ToneDelivery(s.MBED.HiTrials);
    timeLow = s.MBED.ToneDelivery(s.MBED.LowTrials);
    timeRew = s.MBED.RewTimes;
    try
        timeLicks = s.MBED.Licks(:,1);
    catch
        disp('No Licks')
        timeLicks = [];
    end

    if length(timeHi) == length(timeRew)
        disp('Number of Tones and Rewards Matched')
    else
        error('Mismatch in Tones and Rewards')
    end

    for i = 1:length(timeHi)
        %first, find the duration of the pre-delivery period
        prePeriod = timeRew(i) - timeHi(i);

        licksPreHi(i) = length(find(timeLicks > timeHi(i) - prePeriod & timeLicks < timeHi(i)));
        licksAntiHi(i) = length(find(timeLicks > timeHi(i) & timeLicks < timeRew(i)));
        licksPostHi(i) = length(find(timeLicks > timeRew(i) & timeLicks < timeRew(i) + 2000));

        licksPreLow(i) = length(find(timeLicks > timeLow(i) - prePeriod & timeLicks < timeLow(i)));
        licksAntiLow(i) = length(find(timeLicks > timeLow(i) & timeLicks < timeLow(i)+prePeriod));
        licksPostLow(i) = length(find(timeLicks > timeLow(i)+prePeriod & timeLicks < timeLow(i)+prePeriod+2000));
    end
    %store
    bigMaster(masterInd:masterInd + length(timeHi) -1,1) = licksPreHi;
    bigMaster(masterInd:masterInd + length(timeHi) -1,2) = licksAntiHi;
    bigMaster(masterInd:masterInd + length(timeHi) -1,3) = licksPostHi;
    bigMaster(masterInd:masterInd + length(timeHi) -1,4) = licksPreLow;
    bigMaster(masterInd:masterInd + length(timeHi) -1,5) = licksAntiLow;
    bigMaster(masterInd:masterInd + length(timeHi) -1,6) = licksPostLow;

    meanLicksPreHi = mean(licksPreHi);
    meanLicksAntiHi = mean(licksAntiHi);
    meanLicksPostHi = mean(licksPostHi);

    meanLicksPreLow = mean(licksPreLow);
    meanLicksAntiLow = mean(licksAntiLow);
    meanLicksPostLow = mean(licksPostLow);
    
    %store
    sumMaster(bigInd,1) = meanLicksPreHi;
    sumMaster(bigInd,2) = meanLicksAntiHi;
    sumMaster(bigInd,3) = meanLicksPostHi;
    sumMaster(bigInd,4) = meanLicksPreLow;
    sumMaster(bigInd,5) = meanLicksAntiLow;
    sumMaster(bigInd,6) = meanLicksPostLow;
    
    %extract traces and store
    
    bigStruct.HiRaster(:,masterInd:masterInd+length(timeHi)-1) = s.PhotoRaster.ToneRaster(:,s.MBED.HiTrials);
    bigStruct.LowRaster(:,masterInd:masterInd+length(timeHi)-1) = s.PhotoRaster.ToneRaster(:,s.MBED.LowTrials);
    bigStruct.RewRaster(:,masterInd:masterInd+length(timeHi)-1) = s.PhotoRaster.RewardRaster;
    
    bigStruct.VelHi(:,masterInd:masterInd+length(timeHi)-1) = s.VelRaster.ToneRaster(:,s.MBED.HiTrials);
    bigStruct.VelLow(:,masterInd:masterInd+length(timeHi)-1) = s.VelRaster.ToneRaster(:,s.MBED.LowTrials);
    bigStruct.VelRew(:,masterInd:masterInd+length(timeHi)-1) = s.VelRaster.RewardRaster;

    %plot out some basic stuff

    figure
    subplot(3,3,1)
    hist(diff(timeITI)/1000)
    title('Histogram of Tone ITIs')

    subplot(3,3,2)
    hist((timeRew - timeHi)/1000);
    title('Histogram of Reward Delay')
    
    subplot(3,3,3)
    title(testNames{bigInd})

    subplot(3,3,4)
    plot(licksPreLow,'r')
    hold on
    plot(licksPreHi)
    title('Licks Pre (r = low b = hi)')

    subplot(3,3,5)
    plot(licksAntiLow,'r')
    hold on
    plot(licksAntiHi)
    title('Licks Anti (r = low b = hi)')

    subplot(3,3,6)
    plot(licksPostLow,'r')
    hold on
    plot(licksPostHi)
    title('Licks Post (r = low b = hi)')
    
    
    masterInd = masterInd+length(timeHi);
end



%plot out some basics

figure
subplot(3,3,1)
plot(sumMaster(:,2) - sumMaster(:,1))
hold on
plot(sumMaster(:,5) - sumMaster(:,4),'r')
title('Anticipatory Licks (corrected for baseline)')
ylabel('Licks')
xlabel('Days of Training')

subplot(3,3,2)
plot(sumMaster(:,2)./sumMaster(:,5))
title('Ratio of Anticipatory Licks')
ylabel('Ratio')
xlabel('Days of Training')

subplot(3,3,3)
plot(sumMaster(:,2) - sumMaster(:,5))
title('Difference in Anticipatory Licks')
ylabel('Lick Difference')
xlabel('Days of Training')

subplot(3,3,4)
plot(bigMaster(:,2))
hold on
plot(bigMaster(:,5),'r')
title('Anticipatory Licks Across Trials')
ylabel('Licks')
xlabel('Trials')

subplot(3,3,5)
hist(bigMaster(:,2),[1:20])
title('Histogram of AntiLicks Hi')

subplot(3,3,6)
hist(bigMaster(:,5),[1:20])
title('Histogram of AntiLicks Low')

subplot(3,3,7)
plot(bigMaster(:,2) - bigMaster(:,1))
hold on
plot(bigMaster(:,5) - bigMaster(:,4),'r')
title('Anticipatory Licks Across Trials')
ylabel('Licks')
xlabel('Trials')





















