function [s] = functionTuningCurveGenerator(toneReps,toneDur,...
    fs,L,paddingL,prePause,postPauseMin,postPauseMax,startF,...
    endF,octFrac,startdB,enddB,dbSteps,TTLDur, maxDB, rampDur);

warningCheck = (postPauseMin/1000 - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

%this generates a vector with the frequencies that will be used
octRange = log2(endF/startF);
freqs = zeros(octRange/octFrac+1,1);
freqs(1) = startF;
for i = 2:octRange/octFrac+1
    freqs (i) = freqs(i-1)*(2^octFrac);
end

%this generates a vector of all decibel steps
dBs = [startdB:-dbSteps:enddB];
amps = ones(length(dBs),1);
for i = 1:length(amps)
    amps(i) = amps(i)*10^-((maxDB-dBs(i))/20);
end

%list of all frequency/dB pairs
fullList = zeros(length(freqs)*length(dBs),3);
counter = 1;

%simply list every combination
for i = 1:length(freqs)
    for j = 1:length(dBs)
        fullList(counter,1) = freqs(i);
        fullList(counter,2) = dBs(j);
        fullList(counter,3) = amps(j);
        counter = counter+1;
    end
end
%computes length of list
listLength = length(fullList);
%check this is correct
if listLength ~= length(freqs)*length(dBs)
    error('ListLength Incorrect')
end

%generate a n x 1 vector that indicates all indices for every tuning trial.
listDesig = zeros(listLength*toneReps,1);
counter = 1;
for i = 1:toneReps
    listDesig(counter:counter + listLength - 1) = randperm(listLength);
    counter = counter + listLength;
end

%creates array for storing everything. Also fills index with information
%from above.
master = zeros(length(freqs)*length(dBs)*toneReps,4);
master(:,1:3) = fullList(listDesig,1:3);


%generates pauses with exponential function

k = 2.5;
p = (1-exp(-k))*rand(length(master),1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 

master(:,4) = x/1000;

%ramp times for onset and offset in seconds
onRampDur = rampDur*fs; 
offRampDur = rampDur*fs;
remainingPoints = L-onRampDur-offRampDur;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;


%this makes the profile for the TTL signal
ttlSig = zeros(paddingL,1);
ttlSig(1:TTLDur*fs/1000) = 1;


%actual code for performing tuning
for i = 1:length(master)
    pause(prePause)
    toneFreq = master(i,1);
    toneAmpl = master(i,3);
    toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
    finalWave = toneWave.*rampProfile;
    finalWave = finalWave*toneAmpl;
    paddedWave = zeros(paddingL,1);
    paddedWave(1:size(finalWave,1)) = finalWave;
    soundVector = [paddedWave,ttlSig];
    sound(soundVector,fs);
    length(master) - i
    pause(master(i,4))
end

soundData = struct;
soundData.Frequencies = master(:,1);
soundData.dBs = master(:,2);
soundData.Amplitudes = master(:,3);
soundData.Delays = master(:,4);
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;
soundData.OnRamp = onRampDur;
soundData.OffRamp = offRampDur;

s=soundData;

end