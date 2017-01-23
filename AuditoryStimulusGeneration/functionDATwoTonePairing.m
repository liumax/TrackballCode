%this is meant to be the function version of the
%soundGeneratorAlternatingToneWithOptoPulsing. 
function [s] = functionDATwoTonePairing(targetFreq,controlFreq,...
    fs,targetDB,controlDB,toneReps,interRepMin,interRepMax,toneDur,laserLag,...
    ttlDur,ttlITI,maxDB,rampDur);

timeBuffer = 0.1;
longTimeBuffer = 0.2;

calcDBtarget = targetDB - maxDB;
calcDBcontrol = controlDB - maxDB;

targetAmpl = 1*10^(calcDBtarget/20);
controlAmpl = 1*10^(calcDBcontrol/20);

onRampDur = rampDur*fs; 
offRampDur = rampDur*fs;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
L = round(fs*toneDur);

rampProfile = ones(L,1);
%generate the overall ramp profile
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;
%% Calculate laser lag and incorporate into duration of sound file

%first, compensate for loss of time due to second TTL pulse
laserLag = laserLag - ttlITI;

if isnumeric(laserLag) && laserLag == 0
    %set length of tone file
    L = round(fs*toneDur);
    paddingL = round(L + round(fs*timeBuffer)); %add buffer
    
    %this makes the profile for the TTL signal
    ttlSig = zeros(paddingL,1);
    ttlSig(1:round(ttlDur*fs)) = 1;
    
    %make laser TTL
    laserTTL = zeros(paddingL,1);
    laserTTL(1:round(ttlDur*fs)) = 1;
    laserTTL(round(fs*ttlITI):round(fs*(ttlITI+ttlDur))) = 1;
    
    %adjust ramp profile
    newRamp = zeros(paddingL,1);
    newRamp(round(ttlDur*fs):round(ttlDur*fs)+length(rampProfile)-1) = rampProfile;
    rampProfile = newRamp;
elseif isnumeric(laserLag) && laserLag > 0
    if laserLag < toneDur
        %set length of tone file
        L = round(fs*toneDur);
        paddingL = round(L + round(fs*timeBuffer)); %add buffer

        %this makes the profile for the TTL signal
        ttlSig = zeros(paddingL,1);
        ttlSig(1:round(ttlDur*fs)) = 1;

        %make laser TTL
        laserTTL = zeros(paddingL,1);
        laserTTL(round(fs*laserLag):round(fs*(laserLag+ttlDur))) = 1;
        laserTTL(round(fs*(laserLag+ttlITI)):round(fs*(laserLag+ttlITI+ttlDur))) = 1;
        %insert tone TTL
        laserTTL(1:round(ttlDur*fs)) = 1;
        
        %reset ramp profile
        newRamp = zeros(paddingL,1);
        newRamp(1:length(rampProfile)) = rampProfile;
        rampProfile = newRamp;

    elseif laserLag >= toneDur
        %set length of tone file
        L = round(fs*laserLag);
        paddingL = round(L + round(fs*longTimeBuffer)); %add 0.2 second buffer

        %this makes the profile for the TTL signal
        ttlSig = zeros(paddingL,1);
        ttlSig(1:round(ttlDur*fs)) = 1;

        %make laser TTL
        laserTTL = zeros(paddingL,1);
        laserTTL(round(fs*laserLag):round(fs*(laserLag+ttlDur))) = 1;
        laserTTL(round(fs*(laserLag+ttlITI)):round(fs*(laserLag+ttlITI+ttlDur))) = 1;
        %insert tone TTL
        laserTTL(1:round(ttlDur*fs)) = 1;
        %reset ramp profile
        newRamp = zeros(paddingL,1);
        newRamp(1:length(rampProfile)) = rampProfile;
        rampProfile = newRamp;
    end
elseif isnumeric(laserLag) && laserLag < 0
    %set length of tone file
    L = round(fs*(toneDur-laserLag));
    paddingL = round(L + round(fs*longTimeBuffer)); %add 0.2 second buffer
    
    %this makes the profile for the TTL signal
    ttlSig = zeros(paddingL,1);
    ttlSig(round(fs*-laserLag):round(fs*(-laserLag+ttlDur))) = 1;
    
    %make laser TTL
    laserTTL = zeros(paddingL,1);
    laserTTL(1:round(ttlDur*fs)) = 1;
    laserTTL(round(fs*ttlITI):round(fs*(ttlITI+ttlDur))) = 1;
    
    %add tone TTL to laser TTL signal
    laserTTL(round(fs*-laserLag):round(fs*(-laserLag+ttlDur))) = 1;
    
    %reset ramp profile
    newRamp = zeros(paddingL,1);
    newRamp(round(fs*-laserLag):round(fs*-laserLag)-1+length(rampProfile)) = rampProfile;
    rampProfile = newRamp;
elseif ~isnumeric(laserLag)
    %set length of tone file
    L = round(fs*toneDur);
    paddingL = round(L + fs*timeBuffer); %add buffer
    
    %this makes the profile for the TTL signal
    ttlSig = zeros(paddingL,1);
    ttlSig(1:round(ttlDur*fs)) = 1;
    
    %make laser TTL
    laserTTL = zeros(paddingL,1);
    laserTTL(1:round(ttlDur*fs)) = 1;
end



%actual code for running behavior!
toneWave = sin(2*pi*(targetFreq/fs)*(1:paddingL))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
targetVector = [paddedWave,laserTTL];

toneWave = sin(2*pi*(controlFreq/fs)*(1:paddingL))';
finalWave = (toneWave.*rampProfile)*controlAmpl;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
controlVector = [paddedWave,ttlSig];

%generates alternating array of ones and zeros
controller = zeros(2*toneReps,1);
controller(1:2:end) = 1;

%generate pauses with exponential distribution
k = 2.5;
p = (1-exp(-k))*rand(length(controller),1);
tau = (interRepMax-interRepMin)/k;
x = (interRepMin + (-log(1-p))*tau); 


for i=1:length(controller)
    if controller(i) == 1
        sound(targetVector,fs);
    elseif controller(i) == 0
        sound(controlVector,fs);
    end
    length(controller)-i
    pause(x(i))
end

freqRecord = controller;
freqRecord(freqRecord == 1) = targetFreq;
freqRecord(freqRecord == 0) = controlFreq;

dbRecord = controller;
dbRecord(dbRecord == 1) = targetDB;
dbRecord(dbRecord == 0) = controlDB;

ampRecord = controller;
ampRecord(ampRecord == 1) = targetAmpl;
ampRecord(ampRecord == 0) = controlAmpl;


soundData = struct;
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;
soundData.ITI = x;
soundData.Frequencies = freqRecord;
soundData.dBs = dbRecord;
soundData.Amplitudes = ampRecord;
soundData.OptoStimDelay = laserLag;
soundData.OnRampDuration = onRampDur;
soundData.OffRampDuration = offRampDur;
soundData.TargetFreq = targetFreq;
soundData.TargetDB = targetDB;
soundData.ControlFreq = controlFreq;
soundData.ControlDB = controlDB;
soundData.LaserTriggerPulseITI = ttlITI;

s = soundData;

end