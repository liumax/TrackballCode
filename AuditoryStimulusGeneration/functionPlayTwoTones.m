%this is meant to be the function version of the
%soundGeneratorAlternatingToneWithOptoPulsing. 
function [s] = functionPlayTwoTones(targetFreq,controlFreq,...
    fs,targetDB,controlDB,toneReps,interRep,toneDur,TTLDur,maxDB);

calcDBtarget = targetDB - maxDB;
calcDBcontrol = controlDB - maxDB;

targetAmpl = 1*10^(calcDBtarget/20);
controlAmpl = 1*10^(calcDBcontrol/20);

onRampDur = 0.1*fs; 
offRampDur = 0.1*fs;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;

L = fs*toneDur;
paddingL = L + fs*0.4;

rampProfile = ones(L,1);
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;
%if statements to calculate times based on different conditions


ttlSig = zeros(paddingL,1);
ttlSig(1:TTLDur*fs/1000) = 1;

interRep = paddingL/fs+interRep;

%actual code for running behavior!
toneWave = sin(2*pi*(targetFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*targetAmpl;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
targetVector = [paddedWave,ttlSig];

toneWave = sin(2*pi*(controlFreq/fs)*(1:L))';
finalWave = (toneWave.*rampProfile)*controlAmpl;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;
controlVector = [paddedWave,ttlSig];

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

dbRecord = controller;
dbRecord(dbRecord == 1) = targetDB;
dbRecord(dbRecord == 0) = controlDB;

ampRecord = controller;
ampRecord(ampRecord == 1) = targetAmpl;
ampRecord(ampRecord == 0) = controlAmpl;

soundData = struct;
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;
soundData.ITI = interRep;
soundData.ToneDuration = toneDur;
soundData.Frequencies = freqRecord;
soundData.dBs = dbRecord;
soundData.Amplitudes = ampRecord;
soundData.OnRamp = onRampDur;
soundData.OffRamp = offRampDur;
soundData.TargetFreq = targetFreq;
soundData.TargetDB = targetDB;
soundData.ControlFreq = controlFreq;
soundData.ControlDB = controlDB;

s = soundData;

end