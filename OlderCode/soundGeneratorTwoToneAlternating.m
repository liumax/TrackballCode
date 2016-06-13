targetFreq = 6727; %target frequency in Hz
controlFreq = 16000;
targetAmpl = 1; %target amplitude as fraction of 1
toneReps = 100; %tone repetitions
interRep = 5; %seconds between tones

toneDur = 5; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs; %number of samples at correct sampling frequency

interRep = toneDur+interRep;

%ramp times for onset and offset in seconds
onRampDur = 0.005; 
offRampDur = 0.005;

controller = zeros(2*toneReps,1);
controller(1:2:end) = 1;

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
targetVector = [finalWave,ttlSig];

toneWave = sin(2*pi*(controlFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
controlVector = [finalWave,ttlSig];


for i=1:length(controller)
    if controller(i) == 1
        sound(targetVector,fs);
    elseif controller(i) == 0
        sound(controlVector,fs);
    end
    pause(interRep)
end

