

%This code is meant to extract further information from two tone pavlov
%basic analysis

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);



testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);

bigMaster = zeros(100,10);
sumMaster = zeros(5,2);
bigStruct = struct;
masterInd = 1;

trialNumStore = zeros(5,2);


for bigInd = 1:length(testNames)
    trialNumStore(bigInd(1),1) = masterInd;
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
    
    bigStruct.Raw470{bigInd} = s.Photo.Photo.x70dF;
    bigStruct.RawTimes{bigInd} = s.Photo.Photo.x70dFTime;
    bigStruct.ToneTimes{bigInd} = s.Photo.MBEDSig;
    bigStruct.RewTimes{bigInd} = interp1(s.MBED.Jitter,s.Photo.Jitter,s.MBED.RewTimes);
    bigStruct.ToneBig{bigInd} = s.MBED.HiTrials;
    if length(s.MBED.Licks > 0)
        bigStruct.Licks{bigInd} = interp1(s.MBED.Jitter,s.Photo.Jitter,s.MBED.Licks);
        bigStruct.Licks{bigInd}(isnan(bigStruct.Licks{bigInd})) = [];
    end
    
    bigStruct.LicksHi(:,bigInd) = s.Licking.ToneHistHi;
    bigStruct.LicksLow(:,bigInd) = s.Licking.ToneHistLow;
    bigStruct.LickTimes(:,bigInd) = s.Licking.Axis;
    

    %plot out some basic stuff

%     figure
%     subplot(3,3,1)
%     hist(diff(timeITI)/1000)
%     title('Histogram of Tone ITIs')
% 
%     subplot(3,3,2)
%     hist((timeRew - timeHi)/1000);
%     title('Histogram of Reward Delay')
%     
%     subplot(3,3,3)
%     title(testNames{bigInd})
% 
%     subplot(3,3,4)
%     plot(licksPreLow,'r')
%     hold on
%     plot(licksPreHi)
%     title('Licks Pre (r = low b = hi)')
% 
%     subplot(3,3,5)
%     plot(licksAntiLow,'r')
%     hold on
%     plot(licksAntiHi)
%     title('Licks Anti (r = low b = hi)')
% 
%     subplot(3,3,6)
%     plot(licksPostLow,'r')
%     hold on
%     plot(licksPostHi)
%     title('Licks Post (r = low b = hi)')
    
    trialNumStore(bigInd(1),2) = masterInd+length(timeHi)-1;
    masterInd = masterInd+length(timeHi);
end


save('Results','bigMaster','bigStruct','trialNumStore','sumMaster')


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
xlim([0 20])

subplot(3,3,6)
hist(bigMaster(:,5),[1:20])
title('Histogram of AntiLicks Low')
xlim([0 20])

subplot(3,3,7)
plot(bigMaster(:,2) - bigMaster(:,1))
hold on
plot(bigMaster(:,5) - bigMaster(:,4),'r')
title('Anticipatory Licks Across Trials')
ylabel('Licks')
xlabel('Trials')

subplot(3,3,8)
hist(bigMaster(:,2) - bigMaster(:,1),[-10:20])
title('Histogram of AntiLicks Hi BaseSub')
xlim([-10 20])

subplot(3,3,9)
hist(bigMaster(:,5) - bigMaster(:,4),[-10:20])
title('Histogram of AntiLicks Low BaseSub')
xlim([-10 20])

%lets take day 9 as a test case, since this has the biggest range in terms
%of some trials in which the no-rew tone produces lick responses, and a
%bunch of times it doesnt.

figure
plot(bigMaster(751:850,2))
hold on
plot(bigMaster(751:850,5),'r')


falseStarts = find(bigMaster(:,5) >= 4);

targetStarts = find(falseStarts >= 751 & falseStarts <= 850);

targetTraces = bigStruct.LowRaster(:,falseStarts(targetStarts));

%find matched number of low licks

noStarts = find(bigMaster(:,5) == 0);
noStarts = noStarts(noStarts >= 751 & noStarts <= 850);
selectionInd = randperm(length(noStarts),length(targetStarts));
noStartTraces = bigStruct.LowRaster(:,noStarts(selectionInd));

meanTarget = mean(targetTraces');
meanNoTarget = mean(noStartTraces');

figure
plot(meanTarget)
hold on
plot(meanNoTarget,'r')

%hmm, looks like a bit of a difference 

%try plotting out while excluding the first trials
meanTarget = mean(targetTraces(:,2:end)');
meanNoTarget = mean(noStartTraces(:,2:end)');
figure
plot(meanTarget)
hold on
plot(meanNoTarget,'r')

%okay, difference went away. LOL.

%hmm...maybe i should work on the peak detection failures instead...

%doing work on the findphotopeaks code now. 

%okay...so made a function that should do this properly. 

%plot all licking behavior over time

figure
hold on
plot(smooth(bigMaster(:,1),21),'b')
plot(smooth(bigMaster(:,2),21),'r')
plot(smooth(bigMaster(:,4),21),'k')
plot(smooth(bigMaster(:,5),21),'g')


%plot out velocities
figure
imagesc(bigStruct.VelHi')
colormap('parula')
figure
imagesc(bigStruct.VelLow')
colormap('parula')

% for i = 1:85
%     averageTraceHi(:,i) = mean(bigStruct.HiRaster(:,(i-1)*10+1:(i*10))');
%     averageTraceLow(:,i) = mean(bigStruct.LowRaster(:,(i-1)*10+1:(i*10))');
%     
% end



 


