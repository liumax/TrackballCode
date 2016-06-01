%This will be master code. The purpose of this code is to do the following:
%First, this code will wait a default 2 minutes for baseline recordings.
%Following this waiting period, the code will perform a tuning curve with
%parameters that are specified in the code. After a brief pause, the code
%will then perform auditory pairing with dopamine terminal/cell body
%stimulation. After this, the code will an extended tuning curve (double
%the number of tone repetitions). This is because of what Schreiner was
%saying about the mouse needing re-exposure to trigger the changes in
%response properties.
%% This is to enter the filename for the recordings that will be produced.
[fname pname] = uiputfile('test1.mat');
periodFinder = strfind(fname,'.');
fileName = fname(1:periodFinder-1);
%%
% general parameters:
fs = 192000; %sampling frequency in Hz
%%
%tuning curve parameters:
tuningReps = 20; %number of repetitions of each tone/amplitude pair
tuningToneDur = 0.1; %tone duration in seconds

tuningL = tuningToneDur*fs; %number of samples at correct sampling frequency
tuningpaddingL = tuningL + fs*0.1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.

prePause = 0.1; %pause in seconds before tone
postPauseMin = 500; %pause in milliseconds after tone
postPauseMax = 1000; %pause in milliseconds after tone

startF = 4000; %starting frequency in Hz
endF = 32000; %ending frequency in Hz
octFrac = 0.5; %fractions of octaves to move

startdB = 100; %starting decibel value
enddB = 40; %lowest decibel value
dbSteps = 20; %resolution of decible steps
%%
%auditory pairing parameters:
targetFreq = 8000; %target frequency in Hz
controlFreq = 16000;
targetAmpl = 1; %target amplitude as fraction of 1
controlAmpl = 1;
baselineToneReps = 20; %tone repetitions for presentations of long tones before pairing
pairingToneReps = 200; %tone repetitions for pairing experiment
interRep = 2; %seconds between tones

pairingToneDur = 1; %tone duration in seconds
optoDelay = 0.6; %delay between tone onset and opto output. Positive means opto follows sound, negative means sound follows opto
optoDur = 1; %duration of all opto pulses.
optoTTL = 0.002; %duration of opto TTL pulse send through audio card.
optoLag = 0.004; %lag due to the double pulse requirement for triggering

%%
%parameters for TTLs that will signal shifts between functions
signalTTLDur = 0.002; %duration of pulses
signalTTLiti = 0.02; %ITI between signal pulse onsets
signalTTLNum = 4; %number of signal pulses
interFunctionPause = 1; %seconds to wait after a function to finish before starting next

%%
%generates signal TTL signal
ttlSig = zeros(fs*0.5,1); %this is to generate TTL signals to the MBED
for i = 1:signalTTLNum
    ttlSig(1+(fs*signalTTLiti*(i-1)):fs*signalTTLDur+(fs*signalTTLiti*(i-1))) = 1;
end
zeroSig = zeros(fs*0.5,1); %this is a zero signal to send to the speaker
signalVector = [zeroSig,ttlSig];

%%
%This part of the code will actually execute the code!!
%This plays the first tuning curve! Saves the file as the FIRST tuning
%curve
functionTuningCurveGenerator(tuningReps,tuningToneDur,...
    fs,tuningL,tuningpaddingL,prePause,postPauseMin,postPauseMax,startF,...
    endF,octFrac,startdB,enddB,dbSteps,strcat(fileName,'First'));

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);

%this will play the stimuli that will soon be paired without any pairing.
%This will help determine if there are any differences in teh response
%properties of cells to longer tones.
functionPlayTwoTones(targetFreq,controlFreq,...
    fs,targetAmpl,controlAmpl,baselineToneReps,...
    interRep,pairingToneDur,strcat(fileName,'First'));

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);

%This should execute the pairing of auditory stimuli with dopamine
%terminal/cellbody stimulation
functionDATwoTonePairing(targetFreq,controlFreq,...
    fs,targetAmpl,controlAmpl,pairingToneReps,interRep,...
    pairingToneDur,optoDelay,...
    optoDur,optoTTL,optoLag,fileName);

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);

%play tones again after pairing.
functionPlayTwoTones(targetFreq,controlFreq,...
    fs,targetAmpl,controlAmpl,baselineToneReps,...
    interRep,pairingToneDur,strcat(fileName,'Second'));

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);

%This plays the second tuning curve! Saves the file as the SECOND tuning
%curve. Plays double the number of repetitions.
functionTuningCurveGenerator(tuningReps*2,tuningToneDur,...
    fs,tuningL,tuningpaddingL,prePause,postPauseMin,postPauseMax,startF,...
    endF,octFrac,startdB,enddB,dbSteps,strcat(fileName,'Second'));


