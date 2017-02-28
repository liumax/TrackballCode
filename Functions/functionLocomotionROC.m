%This function is meant to parse out rotary encoder based locomotion data,
%compare it to the firing rate of a single unit, and then determine the AUC
%for that unit. It also generates a shuffled set of spikes that serve to
%generate a baseline distribution to compare the AUC value against. 

%Inputs: 
%s: structured array that already has rotary data processed, as well as
%individual units.
%targetName: designation name for the target unit for analysis (ex:
%'nt1cluster4')

%Outputs:
%s: just outputs things into s.



function [s] = functionLocomotionROC(s,targetName);

%important parameters
threshVel = 1; %threshold for calculating a locomotion.
smoothTime = 0.5; %time in seconds for smoothing
rateInc = 20; %increments for going through the firing rates to use as criterion!
velDownSamp = 5; %factor to downsample velocity. Baseline is at 10ms bins.
shuffleReps = 1000; %repetitions of spike shuffling

%first, generate baseline information based on shuffled spike times.

%create store for AUC values
aucStore = zeros(shuffleReps,1);

%get spike times and differences in time
spikeTimes = s.(targetName).SpikeTimes;
spikeTimeDiff = diff(spikeTimes);

%generate smoothed velocity
smoothVel = smooth(s.RotaryData.Velocity(:,2),round(smoothTime/mean(diff(s.RotaryData.Velocity(:,1)))));
smoothVel = downsample(smoothVel,velDownSamp);


%calculate locomotion times:
%find all "locomotor" times
locomotionInd = double(smoothVel>threshVel);
trueLocs = length(find(locomotionInd == 1));
truePause = length(find(locomotionInd == 0));

for shuffInd = 1:shuffleReps

    %make random index of length equal to spikeTimeDiff
    shuffleTimes = randperm(length(spikeTimeDiff));
    shuffleTimes = spikeTimeDiff(shuffleTimes);
    
    shuffSpikeTimes = zeros(length(spikeTimes),1);
    shuffSpikeTimes(1) = spikeTimes(1)+rand;
    shuffSpikeTimes(2:end) = cumsum(shuffleTimes)+shuffSpikeTimes(1);

    fireRate = hist(shuffSpikeTimes,[s.RotaryData.Velocity(1,1):mean(diff(s.RotaryData.Velocity(:,1))):s.RotaryData.Velocity(end,1)])/mean(diff(s.RotaryData.Velocity(:,1)));
    fireRate = downsample(fireRate,velDownSamp);
    %remove first and last bins to eliminate hist errors.
    fireRate(1) = 0;
    fireRate(end) = 0;
    %smooth over the same period

    smoothRate = smooth(fireRate,round(smoothTime/mean(diff(s.RotaryData.Velocity(:,1)))));
    %find min and max rates!
    minRate = (min(smoothRate));
    maxRate = (max(smoothRate));

    %rates at which I will threshold as classifier
    rateRange = minRate:(maxRate-minRate)/rateInc:maxRate;

    rocStore = zeros(4,rateInc);

    for i = 1:rateInc
        %need to conver i to the targeted rate
        threshRate = rateRange(i);
        %find all points at which we classify as locomotion
        classifyInd = (double(smoothRate>=threshRate)+1)*2;
        %compare by adding to locomotionInd
        testInd = classifyInd + locomotionInd;
        %find points with locomotion and classifier(4+1 = 5)
        rocStore(i,1) = length(find(testInd == 5));
        %find points with locomotion but no classifier (2+1 = 3)
        rocStore(i,2) = length(find(testInd == 3));
        %find points with no locomotion but classifier (4 + 0 = 4)
        rocStore(i,3) = length(find(testInd == 4));
        %find points with no locomotion and no classifier (2 + 0 = 2)
        rocStore(i,4) = length(find(testInd == 2));
    end

    %calculate AUC using trapz
    falsePos = rocStore(:,3)/truePause;
    truePos = rocStore(:,1)/trueLocs;
    %reorder in order of false positive from 0 to 1
    [B,I] = sort(falsePos);
    falsePos = B;
    truePos = truePos(I);
    %eliminate duplicate values. 
    [C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
    truePos = truePos(ia);
    falsePos = C;
    %calculate estimate of area under curve. 
    aucStore(shuffInd) = trapz(falsePos,truePos);
    
    if rem(shuffInd,100) == 0;
        disp(strcat('Shuffled Trial:',num2str(shuffInd)))
    end
end

%now calculate actual ROC!


fireRate = hist(spikeTimes,[s.RotaryData.Velocity(1,1):mean(diff(s.RotaryData.Velocity(:,1))):s.RotaryData.Velocity(end,1)])/mean(diff(s.RotaryData.Velocity(:,1)));
fireRate = downsample(fireRate,velDownSamp);
%remove first and last bins to eliminate hist errors.
fireRate(1) = 0;
fireRate(end) = 0;
%smooth over the same period

smoothRate = smooth(fireRate,round(smoothTime/mean(diff(s.RotaryData.Velocity(:,1)))));
%find min and max rates!
minRate = (min(smoothRate));
maxRate = (max(smoothRate));

%rates at which I will threshold as classifier
rateRange = minRate:(maxRate-minRate)/rateInc:maxRate;

rocStore = zeros(4,rateInc);

for i = 1:rateInc
    %need to conver i to the targeted rate
    threshRate = rateRange(i);
    %find all points at which we classify as locomotion
    classifyInd = (double(smoothRate>=threshRate)+1)*2;
    %compare by adding to locomotionInd
    testInd = classifyInd + locomotionInd;
    %find points with locomotion and classifier(4+1 = 5)
    rocStore(i,1) = length(find(testInd == 5));
    %find points with locomotion but no classifier (2+1 = 3)
    rocStore(i,2) = length(find(testInd == 3));
    %find points with no locomotion but classifier (4 + 0 = 4)
    rocStore(i,3) = length(find(testInd == 4));
    %find points with no locomotion and no classifier (2 + 0 = 2)
    rocStore(i,4) = length(find(testInd == 2));
end

%calculate AUC using trapz
falsePos = rocStore(:,3)/truePause;
truePos = rocStore(:,1)/trueLocs;
%reorder in order of false positive from 0 to 1
[B,I] = sort(falsePos);
falsePos = B;
truePos = truePos(I);
%eliminate duplicate values. 
[C,ia,ic] = unique(falsePos,'rows'); %MUST BE ROWS
truePos = truePos(ia);
falsePos = C;
%calculate estimate of area under curve. 
trueAUC= trapz(falsePos,truePos);

s.(targetName).TrueAUC = trueAUC;
s.(targetName).ShuffleAUC = aucStore;

end