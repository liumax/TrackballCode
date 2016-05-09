targetFreq = 8000; %target frequency in Hz
targetAmpl = 1; %target amplitude as fraction of 1
toneReps = 1000; %tone repetitions
interRep = 5; %seconds between tones

toneDur = 0.5; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
optoDelay = 0.4; %delay between tone onset and opto output. Positive means opto follows sound, negative means sound follows opto
optoDur = 1; %duration of all opto pulses.
optoLag = 0.004; %lag due to the double pulse requirement for triggering

onRampDur = 0.005*fs; 
offRampDur = 0.005*fs;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;


if optoDelay < 0
    if optoDur > toneDur-optoDelay
        L = fs*(optoDur);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1+optoLag)) = 0;
        rampProfile(fs*(optoDelay*-1+optoLag):fs*(optoDelay*-1+optoLag)+onRampDur+1) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs/1000) = 1;
        ttlSig(4*fs/1000:5*fs/1000) = 1;
        %want single pulse to mark tone start
        ttlSig((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-0.001)*-1*fs) = 1;
    elseif optoDur <= toneDur - optoDelay
        L = fs*(toneDur-optoDelay);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1)) = 0;
        rampProfile(fs*(optoDelay*-1):fs*(optoDelay*-1)+onRampDur+1) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs/1000) = 1;
        ttlSig(4*fs/1000:5*fs/1000) = 1;
        %want single pulse to mark tone start
        ttlSig((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-0.001)*-1*fs) = 1;
    end
elseif optoDelay > 0
    if toneDur > optoDelay + optoDur
        L = fs*toneDur;
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs/1000) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+0.001)*fs) = 1;
        ttlSig((optoDelay-optoLag+0.004)*fs:(optoDelay-optoLag+0.005)*fs) = 1;
    elseif toneDur < optoDelay + optoDur
        L = fs*(optoDelay + optoDur);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs/1000) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+0.001)*fs) = 1;
        ttlSig((optoDelay-optoLag+0.004)*fs:(optoDelay-optoLag+0.005)*fs) = 1;
    end
elseif optoDelay == 0
    if optoDur + optoLag >= toneDur
        L = fs*(optoDur + optoLag);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoLag)) = 0;
        rampProfile(fs*(optoLag):fs*(optoLag)+onRampDur+1) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
    elseif optoDur + optoLag < toneDur
        L = fs*(toneDur);
        rampProfile = ones(L,1);
    end
end

interRep = toneDur+interRep;
%ramp times for onset and offset in seconds
onRampDur = 0.005*fs; 
offRampDur = 0.005*fs;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;


%this makes the profile for the TTL signal
ttlSig = zeros(L,1);
ttlSig(1:fs/1000) = 1;
ttlSig(4*fs/1000:5*fs/1000) = 1;

%actual code for running behavior!
toneWave = sin(2*pi*(targetFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
soundVector = [finalWave,ttlSig];

for i=1:toneReps
    sound(soundVector,fs);
    pause(interRep)
end

