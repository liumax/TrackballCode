[fname pname] = uiputfile('test1.mat');

toneReps = 1; %number of repetitions of each tone/amplitude pair
toneDur = 0.05; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency

prePause = 0.1; %pause in seconds before tone
postPause = 0.2; %pause in seconds after tone

warningCheck = (postPause - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

startF = 4000; %starting frequency in Hz
endF = 64000; %ending frequency in Hz
octFrac = 0.5; %fractions of octaves to move

startdB = 80; %starting decibel value
enddB = 30; %lowest decibel value
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
amps = zeros (length(dBs),1);
amps(1) = 1;
for i = 2:length(amps)
    amps(i) = (amps(i-1)/sqrt(10));
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

%ramp times for onset and offset in seconds
onRampDur = 0.005; 
offRampDur = 0.005;

%this code generates linear ramps for onset and offset. this reduces issues
%with sharp white noises
rampProfile = ones(L,1);
rampProfile(1:(onRampDur*fs)) = [0:1/(onRampDur*fs):1-1/(onRampDur*fs)];
rampProfile(end-(onRampDur*fs):end) = [1:-1/(onRampDur*fs):0];


%this makes the profile for the TTL signal
ttlSig = zeros(L,1);
ttlSig(1:fs/1000) = 1;

%actual code for performing tuning
for i = 1:length(master)
    pause(prePause)
    toneFreq = master(i,1);
    toneAmpl = master(i,3);
    toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
    finalWave = toneWave.*rampProfile;
    finalWave = finalWave*toneAmpl;
    soundVector = [finalWave,ttlSig];
    sound(soundVector,fs);
    pause(postPause)
end

soundData = struct;
soundData.Frequencies = master(:,1);
soundData.dBs = master(:,2);
soundData.Amplitudes = master(:,3);
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;



save(fullfile(pname,fname),'soundData');
