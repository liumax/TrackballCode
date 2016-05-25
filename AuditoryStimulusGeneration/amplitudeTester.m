targetFreq = 80; %target frequency in Hz
targetAmpl = 1; %target amplitude as fraction of 1
toneReps = 1000; %tone repetitions
interRep = 1; %seconds between tones

toneDur =0.1; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency
totalL = L + fs*0.1;

interRep = toneDur+interRep;
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
ttlSig = zeros(totalL,1);
ttlSig(1:fs/1000) = 1;
% ttlSig(4*fs/1000:5*fs/1000) = 1;

%actual code for running behavior!
toneWave = sin(2*pi*(targetFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
bufferWave = zeros(totalL,1);
bufferWave(1:size(finalWave,1)) = finalWave;
soundVector = [bufferWave,ttlSig];

for i=1:toneReps
    tic
    sound(soundVector,fs);
    toc
    pause(interRep)
end

