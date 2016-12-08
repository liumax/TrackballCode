
%pull dataset
% %first tuning dataset
% tester = s.nt8cluster5.Tuning1Analysis.FreqDBRasters;
% %extract rasters specific to first frequency
% test1 = tester{1};
% [C,ia,ic] = unique(test1(:,2));
% test1(:,3) = [];
% test1(:,2) = ic;
% spikes = test1;
% rasterWindow = [-0.2 0.3];
% numTrials = length(C);
% histBin = 0.01;


%first pairing dataset.
tester = s.nt8cluster5.Pairing2Analysis.AllRasters;
spikes = tester;
[C,ia,ic] = unique(spikes(:,2));
spikes(:,2) = ic;
rasterWindow = [-1 2];
numTrials = length(C);
histBin = 0.02;




stepSize = 0.001; %want to increment by 1ms

numBins = (rasterWindow(2) - rasterWindow(1))/histBin;

timeSteps = [rasterWindow(1)+stepSize:stepSize:rasterWindow(2)]; %generate vector of all mini-steps
timeStepVector = [rasterWindow(1) + stepSize/2:stepSize:rasterWindow(2)];
timeDesig = ceil(timeSteps/histBin); %rounds down to generate the index for pulling from histogram
%need to build in checker if raster goes into negative time
if length(find(timeDesig<0)) > 0
    firstTime = timeDesig(1);
    timeDesig = timeDesig - (timeDesig(1) - 1);
end

% spikeEst1 = zeros(length(genHist),1);
eps = 5;
numReps = 100;
spikeEst = zeros(numBins,numTrials,numReps);
convergeCounter = 0;
nReps = 1;

while convergeCounter <2; 
    %determine if this is an even or odd iteration
    evenOdd = mod(nReps,2);
    %if odd, go forward, if even, go backwards
    if evenOdd == 1;
        trialStart = 1;
        trialInc = 1;
        trialEnd = numTrials;
    elseif evenOdd == 0;
        trialStart = numTrials;
        trialInc = -1;
        trialEnd = 1;
    end
    %next, go through trials
    for nTrial = trialStart:trialInc:trialEnd
        %pull spikes for the targeted trial (time points relative to cue)
        spikeTimes = spikes(spikes(:,2) == nTrial,1);
        %pull a histogram to get all spike times (assuming one spike per
        %1ms bin)
        [spikeCounts ~] = hist(spikeTimes,timeStepVector);
        spikeCounts(spikeCounts > 0) = 1; %make all things equal to one
        %size check
        if length(spikeCounts) ~= length(timeSteps)
            error('Spike Counts dont equal time steps')
        end
        %perform actual adaptive estimation
        for k = 1:length(timeSteps)
            spikeEst(timeDesig(k),nTrial,nReps) = spikeEst(timeDesig(k),nTrial,nReps) + eps*(spikeCounts(k) - spikeEst(timeDesig(k),nTrial,nReps)*stepSize);
        end
        %if not an edge case, then copy estimated spikes to next trial
        if evenOdd == 1 & nTrial < trialEnd | evenOdd == 0 & nTrial > trialEnd
            spikeEst(:,nTrial + trialInc,nReps) = spikeEst(:,nTrial,nReps);
        end
    end
    
    %check for convergence
    if nReps <= 2
        if nReps < numReps
            spikeEst(:,:,nReps + 1) = spikeEst(:,:,nReps);
        end
        disp(strcat('Repetition Number ',num2str(nReps)))
        nReps = nReps + 1; %update repetition number
    elseif nReps > 2
        %find result from 2 iterations ago (this would mean they are going
        %in the same direction. Calculate difference as percentage of
        %values from latest run
        skipDiff = (spikeEst(:,:,nReps - 2) - spikeEst(:,:,nReps))./spikeEst(:,:,nReps);
        diffMax = max(max(abs(skipDiff)));
        if diffMax < 0.05
            convergeCounter = convergeCounter + 1;
            if convergeCounter == 1
                disp(strcat('First Convergence Achieved in ',num2str(nReps),' Repetitions'))
                if nReps < numReps
                    spikeEst(:,:,nReps + 1) = spikeEst(:,:,nReps);
                end
                disp(strcat('Repetition Number',num2str(nReps)))
                nReps = nReps + 1; %update repetition number
            elseif convergeCounter == 2
                disp(strcat('Second Convergence Achieved in ',num2str(nReps),' Repetitions'))
            end
        else
            if nReps < numReps
                spikeEst(:,:,nReps + 1) = spikeEst(:,:,nReps);
            end
            disp(strcat('Repetition Number',num2str(nReps)))
            nReps = nReps + 1; %update repetition number
        end
    end
end