%This will be master code. The purpose of this code is to do the following:
%First, this code will wait a default 2 minutes for baseline recordings.
%Following this waiting period, the code will perform a tuning curve with
%parameters that are specified in the code. After a brief pause, the code
%will actually play the stimuli meant for pairing with dopamine
%stimulation, but without actual laser light. The code then performs the
%auditory pairing with dopamine terminal/cell body stimulation. Following
%this, the code repeats the stimuli without laser, and then plays an
%extended tuning curve(double the number of tone repetitions). This is 
%because of what Schreiner was saying about the mouse needing re-exposure 
%to trigger the changes in response properties.
%% This is to enter the filename for the recordings that will be produced.
[fname pname] = uiputfile('OneTonePairingExtended.mat');
periodFinder = strfind(fname,'.');
fileName = fname(1:periodFinder-1);
%%
% general parameters:
fs = 192000; %sampling frequency in Hz
firstWait = 120; %first waiting period in seconds. 
interFunctionPause = 10; %seconds to wait after a function to finish before starting next
postPairPause = 300; %wait after pairing.
TTLDur = 2; %duration of TTL pulses is ms to be sent via the sound card.
maxDB = 100; %maximum DB that the system is set to. 
rampDur = 0.005; %duration of ramp for tone, in seconds!!!
%%
%tuning curve parameters:
tuningReps = 100; %number of repetitions of each tone/amplitude pair
secondTuningFactor = 5; %factor of extra repetitions for the final tuning after pairing. 
toneDur = 0.1; %tone duration in seconds

tuningL = toneDur*fs; %number of samples at correct sampling frequency
tuningpaddingL = tuningL + fs*0.1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.

prePause = 0.1; %pause in seconds before tone
postPauseMin = 800; %pause in milliseconds after tone
postPauseMax = 1400; %pause in milliseconds after tone

startF = 4000; %starting frequency in Hz
endF = 32000; %ending frequency in Hz
octFrac = 1; %fractions of octaves to move

startdB = 100; %starting decibel value
enddB = 100; %lowest decibel value
dbSteps = 20; %resolution of decible steps
%%
%auditory pairing parameters:
targetFreq = 16000; %target frequency in Hz
targetDB = 100; %target DBs. 100 is max.
pairingToneReps = 200; %tone repetitions for pairing experiment
interRep = 3; %seconds between tones

optoDelay = 0.3; %delay between tone onset and opto output. Positive means opto follows sound, negative means sound follows opto
optoDur = 1; %duration of all opto pulses, in seconds. THIS IS SHITTY HARD CODED VALUE SHOULD EVENTUALLY CHANGE
optoTTL = 0.002; %duration of opto TTL pulse send through audio card.
optoLag = 0.004; %lag due to the double pulse requirement for triggering

%% Confirm frequencies okay
if targetFreq < 4000
    error('Target Frequency Below 4 kHz Cutoff')
elseif startF < 4000
    error('Starting Frequency Below 4 kHz Cutoff')
end

%%
%parameters for TTLs that will signal shifts between functions
signalTTLDur = TTLDur/1000;
signalTTLiti = 0.02; %ITI between signal pulse onsets
signalTTLNum = 4; %number of signal pulses


%%test
%generates signal TTL signal
ttlSig = zeros(fs*0.5,1); %this is to generate TTL signals to the MBED
for i = 1:signalTTLNum
    ttlSig(1+(fs*signalTTLiti*(i-1)):fs*signalTTLDur+(fs*signalTTLiti*(i-1))) = 1;
end
zeroSig = zeros(fs*0.5,1); %this is a zero signal to send to the speaker
signalVector = [zeroSig,ttlSig];

%%
%generates a structured array for storage of all data
fullData = struct;

%% Generate names for all time periods!
names = cell(3,1);
names{1} = 'Tuning1';
names{2} = 'Tuning2';
names{3} = 'Pairing';


%%
%This part of the code will actually execute the code!!
%This pauses the program for the set waiting period to get baseline
%information
disp('Begin Wait Period')
pause(firstWait);
disp('End Wait Period')
%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);
pause(interFunctionPause);

%This plays the first tuning curve! Saves the file as the FIRST tuning
%curve
[s] = functionTuningCurveGenerator(tuningReps,toneDur,...
    fs,tuningL,tuningpaddingL,prePause,postPauseMin,postPauseMax,startF,...
    endF,octFrac,startdB,enddB,dbSteps,TTLDur,maxDB,rampDur);

fullData.(names{1})=s;
s = [];

disp('First Tuning Complete')

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);
pause(interFunctionPause);

%This should execute the pairing of auditory stimuli with dopamine
%terminal/cellbody stimulation
[s] = functionDAOneTonePairing(targetFreq,...
    fs,targetDB,pairingToneReps,interRep,toneDur,optoDelay,...
    optoDur,optoTTL,optoLag,maxDB, rampDur);

fullData.(names{3})=s;
s = [];

disp('Pairing Complete')

%pauses and sends out signal TTLs
pause(postPairPause);

disp('Post Pairing Pause Ended')

sound(signalVector,fs);
pause(interFunctionPause);

%This plays the second tuning curve! Saves the file as the SECOND tuning
%curve. Plays double the number of repetitions.
[s] = functionTuningCurveGenerator(tuningReps*secondTuningFactor,toneDur,...
    fs,tuningL,tuningpaddingL,prePause,postPauseMin,postPauseMax,startF,...
    endF,octFrac,startdB,enddB,dbSteps,TTLDur,maxDB,rampDur);

fullData.(names{2})=s;
s = [];

disp('Second Tuning Complete')

fullData.AllTTLDuration = TTLDur;
fullData.BaselineDuration = firstWait;
fullData.InterfunctionPausing = interFunctionPause;
fullData.SamplingFrequency = fs;

fullData.DividerTTLNumber = signalTTLNum;
fullData.DividerTTLiti = signalTTLiti;

%to calculate number of repetitions in tuning curve
octRange = log2(endF/startF);
freqs = (octRange/octFrac+1);
dBs = size([startdB:-dbSteps:enddB],2);
totalToneReps = freqs*dBs*tuningReps;

fullData.TuningRepetitions = totalToneReps;
fullData.PairingRepetitions = pairingToneReps;
fullData.secondTuningFactor = secondTuningFactor;
fullData.Names = names;
fullData.RampDuration = rampDur;


save(fullfile(pname,fname),'fullData');

