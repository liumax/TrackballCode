%this is meant to be the function version of the
%soundGeneratorAlternatingToneWithOptoPulsing. 
function [] = functionDATwoTonePairing(targetFreq,controlFreq,...
    fs,targetAmpl,controlAmpl,toneReps,interRep,toneDur,optoDelay,...
    optoDur,optoTTL,optoLag,fileName);

pname = pwd;
fname = strcat(fileName,'AuditoryPairing','.mat');

targetFreq = 8000; %target frequency in Hz
controlFreq = 16000;
targetAmpl = 1; %target amplitude as fraction of 1
controlAmpl = 1;
toneReps = 200; %tone repetitions
interRep = 2; %seconds between tones

toneDur = 1; %tone duration in seconds
fs = 192000; %sampling frequency in Hz
optoDelay = 0.6; %delay between tone onset and opto output. Positive means opto follows sound, negative means sound follows opto
optoDur = 1; %duration of all opto pulses.
optoTTL = 0.001; %duration of opto TTL pulse send through audio card.
optoLag = 0.004; %lag due to the double pulse requirement for triggering

onRampDur = 0.1*fs; 
offRampDur = 0.1*fs;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
%if statements to calculate times based on different conditions
if optoDelay < 0 %if opto leads tone
    if optoDur > toneDur-optoDelay %if opto stim exceeds duration of tone and lag period
        L = fs*(optoDur);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1+optoLag)) = 0;
        rampProfile(fs*(optoDelay*-1+optoLag):fs*(optoDelay*-1+optoLag)+onRampDur) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        controlTTL = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*optoTTL) = 1;
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
        %want single pulse to mark tone start
        ttlSig((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-optoTTL)*-1*fs) = 1;
        controlTTL((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-optoTTL)*-1*fs) = 1;
    elseif optoDur <= toneDur - optoDelay %if opto stim is shorter than duration of tone and lag period
        L = fs*(toneDur-optoDelay);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoDelay*-1)) = 0;
        rampProfile(fs*(optoDelay*-1):fs*(optoDelay*-1)+onRampDur) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        controlTTL = zeros(L,1);
        %want double pulse at beginning to trigger laser.
        ttlSig(1:fs*optoTTL) = 1;
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
        %want single pulse to mark tone start
        ttlSig((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-optoTTL)*-1*fs) = 1;
        controlTTL((optoDelay-optoLag)*-1*fs:(optoDelay-optoLag-optoTTL)*-1*fs) = 1;
    end
elseif optoDelay > 0 %in cases where opto follows the tone
    if toneDur > optoDelay + optoDur %if the tone is longer than the opto stim and delay period
        L = fs*toneDur;
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        controlTTL = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        controlTTL(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+optoTTL)*fs) = 1;
        ttlSig((optoDelay)*fs:(optoDelay+optoTTL)*fs) = 1;
    elseif toneDur <= optoDelay + optoDur %if the tone is shorter than the opto stim and delay period
        L = fs*(optoDelay + optoDur);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(fs*toneDur-offRampDur:fs*toneDur) = offRampProfile;
        rampProfile(fs*toneDur:end) = 0;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        controlTTL = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        controlTTL(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig((optoDelay-optoLag)*fs:(optoDelay-optoLag+optoTTL)*fs) = 1;
        ttlSig((optoDelay)*fs:(optoDelay+optoTTL)*fs) = 1;
    end
elseif optoDelay == 0 %if opto stim is coincident with tone
    if optoDur + optoLag >= toneDur %if opto stim and lag is longer than duration of tone
        L = fs*(optoDur + optoLag);
        %GENERATES RAMP AT CORRECT TIME
        rampProfile = ones(L,1);
        rampProfile(1:fs*(optoLag)) = 0;
        rampProfile(fs*(optoLag):fs*(optoLag)+onRampDur) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        controlTTL = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        controlTTL(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
    elseif optoDur + optoLag < toneDur %if tone duration exceeds opto stim
        L = fs*(toneDur);
        rampProfile = ones(L,1);
        rampProfile(1:onRampDur+1) = onRampProfile;
        rampProfile(end-offRampDur:end) = offRampProfile;
        %generates TTL signature. 
        ttlSig = zeros(L,1);
        controlTTL = zeros(L,1);
        %want single pulse at beginning to mark tone
        ttlSig(1:fs*optoTTL) = 1;
        controlTTL(1:fs*optoTTL) = 1;
        %double pulse to trigger laser
        ttlSig(optoLag*fs:(optoLag + optoTTL)*fs) = 1;
    end
end

paddingL = L + fs*0.4;

interRep = paddingL/fs+interRep;

%actual code for running behavior!
toneWave = sin(2*pi*(targetFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
paddedTTL = zeros(paddingL,1);
paddedTTL(1:size(ttlSig,1)) = ttlSig;
targetVector = [paddedWave,paddedTTL];

toneWave = sin(2*pi*(controlFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*controlAmpl;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
paddedTTL = zeros(paddingL,1);
paddedTTL(1:size(controlTTL,1)) = controlTTL;
controlVector = [paddedWave,paddedTTL];

%generates alternating array of ones and zeros
controller = zeros(2*toneReps,1);
controller(1:2:end) = 1;

for i=1:length(controller)
    if controller(i) == 1
        sound(targetVector,fs);
    elseif controller(i) == 0
        sound(controlVector,fs);
    end
    length(controller)-i
    pause(interRep)
end

freqRecord = controller;
freqRecord(freqRecord == 1) = targetFreq;
freqRecord(freqRecord == 0) = controlFreq;

soundData = struct;
soundData.PairedTone.Frequency = targetFreq;
soundData.ControlTone.Frequency= controlFreq;
soundData.PairedTone.Amplitude = targetAmpl;
soundData.ControlTone.Amplitude= controlAmpl;
soundData.Repetitions = toneReps;
soundData.ITI = interRep;
soundData.ToneDuration = toneDur;
soundData.OptoStimDelay = optoDelay;
soundData.ToneOrder = freqRecord;

save(fullfile(pname,fname),'soundData');

end