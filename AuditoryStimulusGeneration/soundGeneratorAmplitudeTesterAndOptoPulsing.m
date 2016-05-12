targetFreq = 8000; %target frequency in Hz
targetAmpl = 1; %target amplitude as fraction of 1
toneReps = 1000; %tone repetitions
interRep = 5; %seconds between tones

mbedDur = 0.001; %duration of pulse to mbed in ms
triggerITI = 0.004; %duration between pulses sent to MBED to trigger laser

toneDur = 0.5; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
optoDelay = -0.1; %delay between tone onset and opto output. Positive means opto follows sound, negative means sound follows opto
optoDur = 1; %duration of all opto pulses.
optoLag = triggerITI; %lag due to the double pulse requirement for triggering

onRampDur = 0.005*fs; 
offRampDur = 0.005*fs;
onRampProfile = (cos((0:1:onRampDur-1)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur-1)/offRampDur*pi)+1)/2;


if optoDelay < 0 %this is if the opto precedes the tone
    if optoDur > toneDur-optoDelay %if optoDuration is greater than tone + delay
        L = round(fs*(optoDur));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1+optoLag)) = 0;
        rampProfile(fs*(optoDelay*-1+optoLag):fs*(optoDelay*-1+optoLag)+onRampDur-1) = onRampProfile;
        rampProfile(fs*(optoDelay*-1+optoLag+toneDur)-offRampDur:fs*(optoDelay*-1+optoLag+toneDur)-1) = offRampProfile;
        rampProfile(fs*(optoDelay*-1+optoLag+toneDur)-1:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*mbedDur) = 1;
        ttlSig(triggerITI*fs:(triggerITI+mbedDur)*fs) = 1;
        %want single pulse to mark tone start
        ttlSig((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-mbedDur)*-1*fs) = 1;
    elseif optoDur <= toneDur - optoDelay %if optoDuration is less than tone + delay
        L = round(fs*(toneDur+(optoDelay*-1)));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1+optoLag)) = 0;
        rampProfile(fs*(optoDelay*-1+optoLag):fs*(optoDelay*-1+optoLag)+onRampDur-1) = onRampProfile;
        rampProfile(end-offRampDur+1:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*mbedDur) = 1;
        ttlSig(triggerITI*fs:(triggerITI+mbedDur)*fs) = 1;
        %want single pulse to mark tone start
        ttlSig(fs*(optoDelay*-1+optoLag):fs*(optoDelay*-1+optoLag+mbedDur)) = 1;
    end
elseif optoDelay > 0 %this is if the opto follows the tone
    if toneDur > optoDelay + optoDur
        L = round(fs*toneDur);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur+1:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs/1000) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+0.001)*fs) = 1;
        ttlSig((optoDelay-optoLag+0.004)*fs:(optoDelay-optoLag+0.005)*fs) = 1;
    elseif toneDur < optoDelay + optoDur
        L = round(fs*(optoDelay + optoDur));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur+1:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs/1000) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+mbedDur)*fs) = 1;
        ttlSig((optoDelay-optoLag+triggerITI)*fs:(optoDelay-optoLag+triggerITI+mbedDur)*fs) = 1;
    end
elseif optoDelay == 0 %This is if opto and tone are supposed to be at the same time
    if optoDur + optoLag >= toneDur
        L = round(fs*(optoDur + optoLag));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoLag)) = 0;
        rampProfile(fs*(optoLag)+1:fs*(optoLag)+onRampDur) = onRampProfile;
        rampProfile(fs*(optoLag)+1+fs*toneDur-offRampDur+1:fs*(optoLag)+1+fs*toneDur) = offRampProfile;
        rampProfile(fs*(optoLag)+1+fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*mbedDur) = 1;
        ttlSig(triggerITI*fs:(triggerITI+mbedDur)*fs) = 1;
        %here, the double pulse also marks the tone!
    elseif optoDur + optoLag < toneDur
        L = round(fs*(toneDur));
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur) = onRampProfile;
        rampProfile(end-offRampDur+1:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*mbedDur) = 1;
        ttlSig(triggerITI*fs:(triggerITI+mbedDur)*fs) = 1;
        %here, the double pulse also marks the tone!
    end
end

interRep = L/fs+interRep;

%actual code for running behavior!
toneWave = sin(2*pi*(targetFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
soundVector = [finalWave,ttlSig];

% for i=1:toneReps
%     sound(soundVector,fs);
%     pause(interRep)
% end

