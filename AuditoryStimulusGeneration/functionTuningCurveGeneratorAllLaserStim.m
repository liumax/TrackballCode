%This is meant to be same code as basic tuning curve, but incorporates
%laser stimulation on every trial.

[fname pname] = uiputfile('test1.mat');

toneReps = 20; %number of repetitions of each tone/amplitude pair
toneDur = 0.1; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency
paddingL = L + fs*0.1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.
ttlTrigDur = 0.002; %time of TTL pulse.
ttlTrigITI = 0.006; %time between TTL pulses for triggering laser.

prePause = 0.1; %pause in seconds before tone
postPauseMin = 600; %pause in milliseconds after tone
postPauseMax = 1000; %pause in milliseconds after tone


warningCheck = (postPauseMin/1000 - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

startF = 4000; %starting frequency in Hz
endF = 32000; %ending frequency in Hz
octFrac = 1; %fractions of octaves to move

startdB = 100; %starting decibel value
enddB = 60; %lowest decibel value
dbSteps = 10; %resolution of decible steps

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
    amps(i) = amps(i)*10^-((maxdB-dBs(i))/20);
end

%list of all frequency/dB pairs
fullList = zeros(length(freqs)*length(dBs),3);
counter = 1;

for i = 1:length(freqs)
    for j = 1:length(dBs)
        fullList(counter:counter+toneReps-1,1) = freqs(i);
        fullList(counter:counter+toneReps-1,2) = dBs(j);
        fullList(counter:counter+toneReps-1,3) = amps(j);
        counter = counter+toneReps;
    end
end

% proxyList = fullList;

fillIndex = zeros(toneReps*length(freqs)*length(dBs),1);
counter = 1;
for i = 1:toneReps*length(dBs)
    fillIndex(counter:counter+length(freqs)-1) = randperm(length(freqs));
    counter = counter + length(freqs);
end

master = zeros(length(freqs)*length(dBs)*toneReps,3);

for i = 1:length(fillIndex)
    master(i,1) = freqs(fillIndex(i));
    holder = find(fullList(:,1) == master(i,1));
    listSize=size(fullList,1);
    holderSize = size(holder,1);
    randIndex = randi(holderSize);
    master(i,2) = fullList(holder(randIndex),2);
    master(i,3) = fullList(holder(randIndex),3);
    fullList(holder(randIndex),:) = [];
end

k = 2.5;
p = (1-exp(-k))*rand(length(master),1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 

master(:,4) = x/1000;

%ramp times for onset and offset in seconds
onRampDur = 0.005*fs; 
offRampDur = 0.005*fs;
remainingPoints = L-onRampDur-offRampDur;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;


%this makes the profile for the TTL signal
ttlSig = zeros(paddingL,1);
ttlSig(1:ttlTrigDur*fs) = 1;
ttlSig(ttlTrigITI*fs:(ttlTrigITI+ttlTrigDur)*fs) = 1;


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



save(fullfile(pname,fname),'soundData');
