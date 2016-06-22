%this is adjusted tuning curve code which should generate a tuning curve
%with alternating modes: half of trials will have laser stim during tone
%delivery, and the other will not. 

[fname pname] = uiputfile('test1.mat');

ttlDur = 0.002; %duration of TTL pulses
ttlITI = 0.006; %ITI for paired pulses.
laserTiming = -0.2; %time in seconds that laser comes relative to tone. negative values mean laser precedes.
prePause = 0.1; %pause in seconds before tone
postPauseMin = 600; %pause in milliseconds after tone
postPauseMax = 1200; %pause in milliseconds after tone

overallReps = 2; %times the whole tuning curve is to be repeated.
toneReps = 20; %number of repetitions of each tone/amplitude pair
toneDur = 0.1; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency
paddingL = L + fs*0.1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.
%laser pairing may require a longer period of sound. 
if laserTiming >= 0 %for laser following tone onset, can just use the same duration of sound data.
    paddingL = L + fs*0.1;
elseif laserTiming < 0
    paddingLPair = L + fs*0.1 + (-laserTiming*fs);
end

warningCheck = (postPauseMin/1000 - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

startF = 4000; %starting frequency in Hz
endF = 32000; %ending frequency in Hz
octFrac = 0.5; %fractions of octaves to move

maxdB = 100; %maximum decibel output
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
%% generates a pseudorandom list of tones and frequencies. 
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
%% calculates intertrial interval.
k = 2.5;
p = (1-exp(-k))*rand(length(master),1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 

master(:,4) = x/1000;

%% repeats all of this to generate a second pseudorandom set of outputs. 
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

pairFillIndex = zeros(toneReps*length(freqs)*length(dBs),1);
counter = 1;
for i = 1:toneReps*length(dBs)
    pairFillIndex(counter:counter+length(freqs)-1) = randperm(length(freqs));
    counter = counter + length(freqs);
end

pairMaster = zeros(length(freqs)*length(dBs)*toneReps,3);

for i = 1:length(pairFillIndex)
    pairMaster(i,1) = freqs(pairFillIndex(i));
    holder = find(fullList(:,1) == pairMaster(i,1));
    listSize=size(fullList,1);
    holderSize = size(holder,1);
    randIndex = randi(holderSize);
    pairMaster(i,2) = fullList(holder(randIndex),2);
    pairMaster(i,3) = fullList(holder(randIndex),3);
    fullList(holder(randIndex),:) = [];
end
%% calculates intertrial interval.
k = 2.5;
p = (1-exp(-k))*rand(length(pairMaster),1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 

pairMaster(:,4) = x/1000;

%% generates ramps for sound to prevent clicks.
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
ttlSig(1:ttlDur*fs) = 1;

ttlSigPairing = zeros(paddingLPair,1);
ttlSigPairing(1:ttlDur*fs) = 1;
ttlSigPairing(ttlITI*fs: (ttlITI+ttlDur)*fs) = 1;
ttlSigPairing(-laserTiming*fs:(-laserTiming+ttlDur)*fs) = 1; %marks tone!


%actual code for performing tuning
for i = 1:length(master)*overallReps;
    pause(prePause)
    
    if mod(i,2) = 0 %calculates remainder of division: this is true for evens
        toneFreq = master(i,1);
        toneAmpl = master(i,3);
        toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
        finalWave = toneWave.*rampProfile;
        finalWave = finalWave*toneAmpl;
        paddedWave = zeros(paddingL,1);
        paddedWave(1:size(finalWave,1)) = finalWave;
        soundVector = [paddedWave,ttlSig];
    else
        toneFreq = pairMaster(i,1);
        toneAmpl = pairMaster(i,3);
        toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
        finalWave = toneWave.*rampProfile;
        finalWave = finalWave*toneAmpl;
        paddedWave = zeros(paddingLPair,1);
        paddedWave(fs*-laserTiming:fs*-laserTiming+size(finalWave,1)-1) = finalWave;
        soundVector = [paddedWave,ttlSigPairing];
    end
    soundVector = [paddedWave,ttlSig];
    sound(soundVector,fs);
    length(master)*overallReps - i
    pause(master(i,4))
end

soundData = struct;
soundData.UnpairedStimuli.Frequencies = master(:,1);
soundData.UnpairedStimuli.dBs = master(:,2);
soundData.UnpairedStimuli.Amplitudes = master(:,3);
soundData.UnpairedStimuli.Delays = master(:,4);
soundData.PairedStimuli.Frequencies = pairMaster(:,1);
soundData.PairedStimuli.dBs = pairMaster(:,2);
soundData.PairedStimuli.Amplitudes = pairMaster(:,3);
soundData.PairedStimuli.Delays = pairMaster(:,4);
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;



save(fullfile(pname,fname),'soundData');
