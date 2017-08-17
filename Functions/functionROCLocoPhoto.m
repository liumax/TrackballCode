% testerCodeApplyingROCFitToPhotometry


function [funcOut] = functionROCLocoPhoto(mbedJitter,tdtJitter,photoTrace,photoTimes,velTrace);


%important parameters
threshVel = 1; %threshold for binary designation for locomotion.
smoothTime = 0.5; %time in seconds for smoothing
rateInc = 20; %increments for going through the firing rates to use as criterion!
velDownSamp = 1; %factor to downsample velocity. Baseline is at 100ms bins.
shuffleReps = 10000; %repetitions of spike shuffling


% velTrace = s.Locomotion.Velocity;
%this has time in MBED timing, need to switch
newVelTime = interp1(mbedJitter/1000,tdtJitter,velTrace(:,1));
velTrace(:,1) = newVelTime;
%remove NANs
velTrace(isnan(velTrace(:,1)),:) = [];
%generate smoothed velocity
smoothVel = smooth(velTrace(:,2),round(smoothTime/mean(diff(velTrace(:,1)))));
% smoothVel = downsample(smoothVel,velDownSamp);



%now we need to process the photometry trace
% photoTrace = s.Photo.Photo.x70dF;
% photoTimes = s.Photo.Photo.x70dFTime;
%match up timing to velocity trace
photoTrace = interp1(photoTimes,photoTrace,velTrace(:,1));

%create store for AUC values
aucStore = zeros(shuffleReps,1);

%calculate locomotion times:
%find all "locomotor" times
locomotionInd = double(smoothVel>threshVel);
trueLocs = length(find(locomotionInd == 1));
truePause = length(find(locomotionInd == 0));

for shuffInd = 1:shuffleReps

    %make random index of length equal to spikeTimeDiff
    shuffleTimes = randperm(length(velTrace(:,1)));
    smoothRate = photoTrace(shuffleTimes);
    
    %find min and max rates!
    minRate = (min(smoothRate));
    maxRate = (max(smoothRate));

    %rates at which I will threshold as classifier
    rateRange = minRate:(maxRate-minRate)/rateInc:maxRate;
    if length(rateRange) == 0
        rocStore = zeros(4,rateInc);
        for i = 1:rateInc
            %need to convert i to the targeted rate
            threshRate = 0;
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
    else

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


smoothRate = photoTrace;
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

funcOut.TrueAUC = trueAUC; %This is the AUC value for the actual spike train
funcOut.ShuffleAUC = aucStore; %these are shuffled spike train AUC values

%works!! XD somehow I got a negative value...

funcOut.TrueAUC = trueAUC; %This is the AUC value for the actual spike train
funcOut.ShuffleAUC = aucStore; %these are shuffled spike train AUC values
funcOut.AUCRange = [prctile(aucStore,0.05),prctile(aucStore,99.95)];

end






