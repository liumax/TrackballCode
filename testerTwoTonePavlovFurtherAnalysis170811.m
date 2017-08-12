

%This code is meant to extract further information from two tone pavlov
%basic analysis

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%first, lets get all of our time points.

%get ITI, tone times, reward delivery times

timeITI = s.MBED.ToneDelivery;
timeHi = s.MBED.ToneDelivery(s.MBED.HiTrials);
timeLow = s.MBED.ToneDelivery(s.MBED.LowTrials);
timeRew = s.MBED.RewTimes;
timeLicks = s.MBED.Licks;

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

meanLicksPreHi = mean(licksPreHi);
meanLicksAntiHi = mean(licksAntiHi);
meanLicksPostHi = mean(licksPostHi);

meanLicksPreLow = mean(licksPreLow);
meanLicksAntiLow = mean(licksAntiLow);
meanLicksPostLow = mean(licksPostLow);


%plot out some basic stuff

figure
subplot(3,3,1)
hist(diff(timeITI)/1000)
title('Histogram of Tone ITIs')

subplot(3,3,2)
hist((timeRew - timeHi)/1000);
title('Histogram of Reward Delay')



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

































