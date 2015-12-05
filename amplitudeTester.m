targetFreq = 8000; %target frequency in Hz
targetAmpl = 1; %target amplitude as fraction of 1


toneDur = 1; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency

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
toneWave = sin(2*pi*(targetFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
soundVector = [finalWave,ttlSig];
sound(soundVector,fs);

