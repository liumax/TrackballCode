
toneReps = 1; %number of repetitions of each tone/amplitude pair
toneDur = 2; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency
paddingL = L + fs*0.1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.

prePause = 0.1; %pause in seconds before tone

pauseTime = 10; %pause in seconds after tone

startF = 4000; %starting frequency in Hz
endF = 64000; %ending frequency in Hz
octFrac = 0.1; %fractions of octaves to move

%this generates a vector with the frequencies that will be used
octRange = log2(endF/startF);
freqs = zeros(octRange/octFrac+1,1);
freqs(1) = startF;
for i = 2:octRange/octFrac+1
    freqs (i) = freqs(i-1)*(2^octFrac);
end

targetDB = 100;
maxDB = 100;
targetAmp = 10^-((maxDB-targetDB)/20);

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

%actual code for performing tuning
for i = 1:length(freqs)
    pause(prePause)
    freqs(i)
    toneFreq = freqs(i);
    toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
    finalWave = toneWave.*rampProfile;
    finalWave = finalWave*targetAmp;
    paddedWave = zeros(paddingL,1);
    paddedWave(1:size(finalWave,1)) = finalWave;
    soundVector = [paddedWave,ttlSig];
    sound(soundVector,fs);
    pause(pauseTime)
end