toneReps = 200;
toneAmp = 0.1;
toneFreq = 4000;
toneDur = 5;
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency
paddingL = L + fs*1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.
% paddingL = L;

pauseTime = 1 + toneDur;

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
ttlSig(1:2*fs/1000) = 1;

toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
finalWave = toneWave.*rampProfile*toneAmp;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
soundVector = [paddedWave,ttlSig];

for i = 1:toneReps
    sound(soundVector,fs);
    pause(pauseTime)
    i
end