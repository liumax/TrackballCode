toneReps = 20; %number of repetitions of each tone/amplitude pair
toneDur = 1; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency

prePause = 0.1; %pause in seconds before tone
postPause = 0.2; %pause in seconds after tone

% laserPrecede = ;
% laserProb = 0; %probability of laser light
% laserOnly = 0; %number of trials with laser only

startF = 4000; %starting frequency in Hz
endF = 64000; %ending frequency in Hz
octFrac = 0.5; %fractions of octaves to move

startdB = 60; %starting decibel value
enddB = 30; %lowest decibel value
dbSteps = 5; %resolution of decible steps

%this generates a vector with the frequencies that will be used
octRange = log2(endF/startF);
freqs = zeros(octRange/octFrac+1,1);
freqs(1) = startF;
for i = 2:octRange/octFrac+1
    freqs (i) = freqs(i-1)*(2^octFrac);
end

%this generates a vector of all decibel steps
dBs = [enddB:dbSteps:startdB];

%list of all frequency/dB pairs
fullList = zeros(length(freqs)*length(dBs),2);
counter = 1;

for i = 1:length(freqs)
    for j = 1:length(dBs)
        fullList(counter:counter+toneReps-1,1) = freqs(i);
        fullList(counter:counter+toneReps-1,2) = dBs(j);
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

master = zeros(length(freqs)*length(dBs)*toneReps,2);

for i = 1:length(fillIndex)
    master(i,1) = freqs(fillIndex(i));
    holder = find(fullList(:,1) == master(i,1));
    listSize=size(fullList,1);
    holderSize = size(holder,1);
    randIndex = randi(holderSize);
    master(i,2) = fullList(holder(randIndex),2);
    fullList(holder(randIndex),:) = [];
end

%%%This code is commented out. previous attempt at making pseudorandom list
%%%of all pairs of amplitudes and frequencies.

% % % master = zeros(length(freqs)*length(dBs)*toneReps,2);
% % % 
% % % %sets first value of master array. This is outside of the inside loop
% % % %because this way I can have different code running this.
% % % counter = 1;
% % % randArray = randperm(length(proxyList));
% % % master(counter,:) = proxyList(randArray(1),:);
% % % proxyList(randArray(1),:) = [];
% % % counter = counter + 1;
% % % 
% % % while ~isempty(proxyList)
% % %     randArray = randperm(size(proxyList,1));
% % %     correctArray = randArray(proxyList(randArray,1)~=master(counter-1,1));
% % %     if isempty(correctArray)
% % %         disp 'FAILURE'
% % %     end
% % %     master(counter,:) = proxyList(correctArray(1),:);
% % %     proxyList(randArray(1),:) = [];
% % %     counter = counter + 1;
% % % end
% % % 
% % % %master array will have frequencies in first column, amplitudes in second
% % % master = zeros(length(freqs)*length(dBs)*toneReps,2);
% % % 
% % % %fills in frequency placeholders
% % % j=1;
% % % for i = 1:length(master)/length(freqs)
% % %     master(j:j+length(freqs)-1,1) = randperm(length(freqs));
% % %     j = j+length(freqs)-1;
% % % end
% % % %replaces with actual frequencies
% % % for i = 1:length(freqs)
% % %     master(master(:,1) == i,1) = freqs(i);
% % % end
% % % %fills in amplitude placeholders
% % % j=1;
% % % for i = 1:length(master)/length(dBs)
% % %     master(j:j+length(dBs)-1,2) = randperm(length(dBs));
% % %     j = j+length(dBs)-1;
% % % end
% % % %replaces with actual amplitudes
% % % for i = 1:length(dBs)
% % %     master(master(:,2) == i,2) = dBs(i);
% % % end
% % % %%%THIS NEEDS ADJUSTING: MUST DO TUNING AND MAP OUT VALUE VS DB

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

%actual code for running behavior!
for i = 1:length(master)
    pause(0.1)
    toneFreq = master(i,1);
    toneAmpl = master(i,2);
    toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
    finalWave = toneWave.*rampProfile;
    soundVector = [finalWave,ttlSig];
    sound(soundVector,fs);
    pause(0.2)
end
