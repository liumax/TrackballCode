[fname pname] = uiputfile('bestFreq.mat');

%This code is generated to do a quick search for best frequency. All tones
%are played at maximum amplitude, with no changes to dB. The whole set of
%tones is repeated in pseudorandom order (finish the full set before doing
%the full set again). Delay has been made into a variable delay to
%hopefully increase salience. 

%tone properties
toneReps = 20; %number of repetitions of each tone/amplitude pair
toneDur = 0.1; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency

%pausing times
prePause = 0.1; %pause in seconds before tone
postPauseMin = 450; %pause in milliseconds after tone
postPauseMax = 1000; %pause in milliseconds after tone

warningCheck = (postPauseMin/1000 - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

%frequency profile
startF = 4000; %starting frequency in Hz
endF = 64000; %ending frequency in Hz
octFrac = 0.5; %fractions of octaves to move

%this generates a vector with the frequencies that will be used
octRange = log2(endF/startF);
freqs = zeros(octRange/octFrac+1,1);
freqs(1) = startF;
for i = 2:octRange/octFrac+1
    freqs (i) = freqs(i-1)*(2^octFrac);
end

%generates pseudorandom indexes of all frequencies
fillIndex = zeros(toneReps*length(freqs),1);
counter = 1;
for i = 1:toneReps
    fillIndex(counter:counter+length(freqs)-1) = randperm(length(freqs));
    counter = counter + length(freqs);
end

%master array to hold major information
master = zeros(length(freqs)*toneReps,2);

for i = 1:length(fillIndex)
    master(i,1) = freqs(fillIndex(i));
end

k = 2.5;
p = (1-exp(-k))*rand(toneReps*length(freqs),1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 
master(:,2) = x/1000;

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
ttlSig(1:fs/100) = 1;

%actual code for performing tuning
for i = 1:length(master)
    pause(prePause)
    toneFreq = master(i,1);
    toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
    finalWave = toneWave.*rampProfile;
    soundVector = [finalWave,ttlSig];
    sound(soundVector,fs);
    length(master) - i
    pause(master(i,2))
end

soundData = struct;
soundData.Frequencies = master(:,1);
soundData.ITI = master(:,2);
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;



save(fullfile(pname,fname),'soundData');
