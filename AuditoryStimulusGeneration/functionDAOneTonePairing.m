%this is meant to be the function version of the
%soundGeneratorAlternatingToneWithOptoPulsing. 
function [s] = functionDAOneTonePairing(targetFreq,...
    fs,targetDB,toneReps,interRep,toneDur,optoDelay,...
    optoDur,optoTTL,optoLag,maxDB, rampDur);

calcDBtarget = targetDB - maxDB;

targetAmpl = 1*10^(calcDBtarget/20);

onRampDur = rampDur*fs; 
offRampDur = rampDur*fs;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
%if statements to calculate times based on different conditions
if optoDelay < 0 %if opto leads tone
    if optoDur > toneDur-optoDelay %if opto stim exceeds duration of tone and lag period
        L = round(fs*(optoDur));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1+optoLag)) = 0;
        rampProfile(fs*(optoDelay*-1+optoLag):fs*(optoDelay*-1+optoLag)+onRampDur) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*optoTTL) = 1;
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
        %want single pulse to mark tone start
        ttlSig((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-optoTTL)*-1*fs) = 1;
    elseif optoDur <= toneDur - optoDelay %if opto stim is shorter than duration of tone and lag period
        L = round(fs*(toneDur-optoDelay));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1)) = 0;
        rampProfile(fs*(optoDelay*-1):fs*(optoDelay*-1)+onRampDur) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*optoTTL) = 1;
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
        %want single pulse to mark tone start
        ttlSig((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-optoTTL)*-1*fs) = 1;
    end
elseif optoDelay > 0 %in cases where opto follows the tone
    if toneDur > optoDelay + optoDur %if the tone is longer than the opto stim and delay period
        L = round(fs*toneDur);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+optoTTL)*fs) = 1;
        ttlSig((optoDelay)*fs:(optoDelay+optoTTL)*fs) = 1;
    elseif toneDur <= optoDelay + optoDur %if the tone is shorter than the opto stim and delay period
        L = round(fs*(optoDelay + optoDur));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+optoTTL)*fs) = 1;
        ttlSig((optoDelay)*fs:(optoDelay+optoTTL)*fs) = 1;
    end
elseif optoDelay == 0 %if opto stim is coincident with tone
    if optoDur + optoLag >= toneDur %if opto stim and lag is longer than duration of tone
        L = round(fs*(optoDur + optoLag));
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoLag)) = 0;
        rampProfile(fs*(optoLag):fs*(optoLag)+onRampDur) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
    elseif optoDur + optoLag < toneDur %if tone duration exceeds opto stim
        L = round(fs*(toneDur));
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
    end
end

paddingL = round(L + fs*0.4);

interRep = paddingL/fs+interRep;

%actual code for running behavior!
toneWave = sin(2*pi*(targetFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
paddedTTL = zeros(paddingL,1);
paddedTTL(1:size(ttlSig,1)) = ttlSig;
targetVector = [paddedWave,paddedTTL];

for i=1:toneReps
    sound(targetVector,fs)
    toneReps - i
    pause(interRep)
end

soundData = struct;
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;
soundData.ITI = interRep;
soundData.ToneDuration = toneDur;
soundData.Frequency = targetFreq;
soundData.dB = targetDB;
soundData.OptoStimDelay = optoDelay;
soundData.OnRampDuration = onRampDur;
soundData.OffRampDuration = offRampDur;
soundData.LaserTriggerPulseITI = optoLag;

s = soundData;

end