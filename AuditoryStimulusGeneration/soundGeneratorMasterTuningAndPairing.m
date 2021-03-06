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
[fname pname] = uiputfile('TuningAndPairing.mat');
periodFinder = strfind(fname,'.');
fileName = fname(1:periodFinder-1);
%% Parameters
%Tuning Curve Parameters
tuningToneDur = 0.1; %tone duration in seconds
tuningReps = 10; %number of repetitions of each tone/amplitude pair

startF = 4000; %starting frequency in Hz
endF = 32000; %ending frequency in Hz
octFrac = 0.5; %fractions of octaves to move

startdB = 100; %starting decibel value
enddB = 100; %lowest decibel value
dbSteps = 20; %resolution of decible steps

secondTuningRatio = 1; %ratio of second tuning/first tuning. >1 means more 
%tones in the second tuning period
prePause = 0.1; %pause in seconds before tone
postPauseMin = 600; %pause in milliseconds after tone
postPauseMax = 1000; %pause in milliseconds after tone

%Tone Presentation Parameters
pairingToneDur = 0.1; %tone duration in seconds
targetFreq = 22627; %target frequency in Hz
controlFreq = 8000;
targetDB = 100; %target DBs. 100 is max.
controlDB = 80;

baselineToneReps = 10; %tone repetitions for presentations of long tones before pairing
pairingToneReps = 10; %tone repetitions for pairing experiment
presPauseMin = 2; %seconds between tones min
presPauseMax = 4;%seconds between tones max

%General Function Parameters:
firstWait = 10; %first waiting period in seconds. 
interFunctionPause = 2; %seconds to wait after a function to finish before starting next

%Laser Pairing Parameters:
optoDelay = 0; %delay between tone onset and opto output. Positive means opto follows sound, negative means sound follows opto
optoTTL = 0.002; %duration of opto TTL pulse send through audio card.
optoLag = 0.004; %lag due to the double pulse requirement for triggering

%% Fixed Parameters
%parameters for TTLs that will signal shifts between functions
TTLDur = 2; %duration of TTL pulses is ms to be sent via the sound card.
signalTTLDur = TTLDur/1000;
signalTTLiti = 0.02; %ITI between signal pulse onsets
signalTTLNum = 4; %number of signal pulses
fs = 192000; %sampling frequency in Hz
maxDB = 100; %maximum DB that the system is set to. 
rampDur = 0.005; %duration of ramp for tone, in seconds!!!


tuningL = tuningToneDur*fs; %number of samples at correct sampling frequency
tuningpaddingL = tuningL + fs*0.1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.

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
names = cell(5,1);
names{1} = 'Tuning1';
names{2} = 'Tuning2';
names{3} = 'Tones1';
names{4} = 'Tones2';
names{5} = 'Pairing';


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
[s] = functionTuningCurveGenerator(tuningReps,tuningToneDur,...
    fs,tuningL,tuningpaddingL,prePause,postPauseMin,postPauseMax,startF,...
    endF,octFrac,startdB,enddB,dbSteps,TTLDur,maxDB,rampDur);

fullData.(names{1})=s;
s = [];

disp('First Tuning Complete')

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);
pause(interFunctionPause);

%this will play the stimuli that will soon be paired without any pairing.
%This will help determine if there are any differences in teh response
%properties of cells to longer tones.
[s] = functionPlayTwoTones(targetFreq,controlFreq,...
    fs,targetDB,controlDB,baselineToneReps,...
    presPauseMin,presPauseMax,pairingToneDur,TTLDur,maxDB,rampDur);

fullData.(names{3})=s;
s = [];

disp('First Tone Presentation Complete')

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);
pause(interFunctionPause);

%This should execute the pairing of auditory stimuli with dopamine
%terminal/cellbody stimulation
[s] = functionDATwoTonePairing(targetFreq,controlFreq,...
    fs,targetDB,controlDB,pairingToneReps,presPauseMin,presPauseMax,...
    pairingToneDur,optoDelay,optoTTL,optoLag,maxDB,rampDur);

fullData.(names{5})=s;
s = [];

disp('Pairing Complete')

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);
pause(interFunctionPause);

%play tones again after pairing.
[s] = functionPlayTwoTones(targetFreq,controlFreq,...
    fs,targetDB,controlDB,baselineToneReps,...
    presPauseMin,presPauseMax,pairingToneDur,TTLDur,maxDB,rampDur);

fullData.(names{4})=s;
s = [];

disp('Second Tone Presentation Complete')

%pauses and sends out signal TTLs
pause(interFunctionPause);
sound(signalVector,fs);
pause(interFunctionPause);

%This plays the second tuning curve! Saves the file as the SECOND tuning
%curve. Plays double the number of repetitions.
[s] = functionTuningCurveGenerator(tuningReps*secondTuningRatio,tuningToneDur,...
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
fullData.SecondTuningRatio = secondTuningRatio;
fullData.PresentationRepetitions = baselineToneReps*2;
fullData.PairingRepetitions = pairingToneReps*2;
fullData.Names = names;
fullData.RampDuration = rampDur;


save(fullfile(pname,fname),'fullData');

