
[fname pname] = uiputfile('toneLaserPair.mat');
%this code is to pair single tones with laser stimulation. Can control gap
%between laser stim and audio delivery. Have control to make series of
%auditory pips per stimulation round. 

%control over tone settings here
toneReps = 150; %number of repetitions of each tone/amplitude pair
totalReps = toneReps;
toneFreq = 8000; %frequency in Hz
toneDur = 0.1; %overall duration in seconds

%control over individual sound pips.
pipNum = 10; %number of tone pips per rep
pipDur = 50; %duration of each pip in msec

%control over laser settings here
% laserDelay = ; %delay between onset of auditory stimulus and laser onset. 
%positive numbers indicate laser follows auditory stimulus, negative
%indicates that laser precedes auditory stimulus.
ttlDur = 0.01; %duration of laser signaling TTL in seconds
ttlGap = 0.005; %delay between laser targeting TTLs
fs = 192000; %sampling frequency in Hz
L = toneDur*fs;

%pausing times!
minPause = 2;
maxPause = 5;

prePause = 0.5; %pause in seconds before tone

%generates random ITIs using exponential function
k = 2.5;
p = (1-exp(-k))*rand(totalReps,1);
tau = (maxPause-minPause)/k;
x = minPause + (-log(1-p))*tau; 


warningCheck = (minPause - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

%generates sound wave.
toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';


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
ttlSig(1:fs*ttlDur) = 1;
ttlSig(fs*ttlDur+fs*ttlGap:fs*ttlDur*2+fs*ttlGap) = 1;

%adds ramps to white noise
finalWave = toneWave.*rampProfile;

%makes two column matrix for sound and TTL output
soundVector = [finalWave,ttlSig];

for i = 1:toneReps
    pause(prePause)
    sound(soundVector,fs);
    disp(toneReps - i)
    pause(x(i))
end

soundData = struct;
soundData.Frequencies = toneFreq;
soundData.ITI = x;
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;
soundData.LaserDelay = ttlGap+ttlDur;

save(fullfile(pname,fname),'soundData');
