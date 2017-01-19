%This code is meant to provide a stimulation script that provides either
%alternating or randomly sequenced sets of two tones, over an extended
%period of time. The paradigm is meant to have three phases: an initial
%baseline phase, in which simply tones are played, followed by a pairing
%phase, in which one tone is paired with optogenetic release of dopamine,
%followed by a final settling phase, in which the tones are played alone



%% Dialog for inputing file name
[fname pname] = uiputfile('ABBA.mat');

%% Parameters

% Tone Parameters
toneDur = 0.1; %tone duration in seconds
onRampDur = 0.005; %onramp duration in seconds Will be cosine
offRampDur = 0.005; %offramp duration in seconds. Will be cosine

toneFreqA = 8000; %toneA will also be the laser targeted tone
toneFreqB = 16000;

toneDBA = 100;
toneDBB = 100;

% Structure Parameters
laserLag = 0.1; %lag of laser. 0 means essentially coincident. -1 means 1 second before, +1 is 1 second after. can also input 'none'
% laserLag = 'none';

baselineReps = 10;
pairingReps = 10;
finishReps = 10;

minITI = 2; %iti in seconds
maxITI = 4;

% orderDes = 'pseudorandom';%can be 'alternating' or 'pseudorandom'
orderDes = 'alternating';


% Fixed Parameters
fs = 192000; %sampling frequency in Hz
prePause = 0.1; %pause in seconds before tone delivery
maxDB = 100;
ttlDur = 0.002;
ttlITI = 0.004; %iti between TTL onsets for laser TTLs
signalTTLiti = 0.02; %ITI between signal pulse onsets
signalTTLNum = 4; %number of signal TTLs
timeBuffer = 0.1; %standard time buffer to pad sound file in seconds.
longTimeBuffer = 0.2; %longer time buffer to pad longer sound files in seconds.

%% Calculate ITIs for tones
numTrials = 2*(baselineReps+pairingReps+finishReps);
k = 2.5;
p = (1-exp(-k))*rand(numTrials,1);
tau = (maxITI-minITI)/k;
x = (minITI + (-log(1-p))*tau); 

master = zeros(numTrials,5);
master(:,1) = x;

%% Safety Checks

%confirm frequencies
if isnumeric(toneFreqA) && toneFreqA < 4000
    error('ToneAFreq is Pure Tone Below 4kHz')
elseif isnumeric(toneFreqB) && toneFreqB < 4000
    error('ToneBFreq is Pure Tone Below 4kHz')
end

%confirm that user wants to continue with protocol duration

totalTime = sum(x)/60; %total time in minutes

prompt = strcat('Protocol Duration is ',num2str(totalTime),'Mins. Proceed (y/n)?');
feedBack = input(prompt,'s');
if strfind(feedBack,'y') 
    disp('Proceeding with Protocol')
elseif strfind(feedBack,'n')
    error('Terminating Protocol on User Feedback')
else
    error('Invalid answer, terminating')
end

%eliminate bad timings for laser timing

if laserLag ~= 0 & abs(laserLag) < ttlITI + ttlDur;
    error('Invalid Laser Lag Time. Terminating')
end

%% Generate order for tone presentations
%order will be 0 for toneA, 1 for toneB
if strfind(orderDes,'alternating')
    orderFile = zeros(numTrials,1);
    orderFile(1:2:end) = 1;
    %check this works out properly
    order1s = length(find(orderFile == 1));
    order0s = length(find(orderFile == 0));
    if order1s ~= order0s
        error('Mismatched Number of Alternating Trials')
    end
elseif strfind(orderDes,'pseudorandom')
    orderFile = randperm(numTrials);
    oddFind = find(mod(orderFile,2) == 1);
    evenFind = find(mod(orderFile,2) == 0);
    %check!
    if length(oddFind) ~= length(evenFind)
        error('Mismatched Number of Pseudorandom Trials')
    end
    orderFile(oddFind) = 0;
    orderFile(evenFind) = 1;
else
    error('Invalid Input for Order')
end

master(:,2) = orderFile;
master(orderFile == 0,3) = toneFreqA;
master(orderFile == 1,3) = toneFreqB;

master(orderFile == 0,4) = toneDBA;
master(orderFile == 1,4) = toneDBB;

%% Calculate some general things for sound presentation
%set default length of tone file
L = round(fs*toneDur);
paddingL = round(fs*timeBuffer); %add buffer

%generate spacer TTLs
spacerTTL = zeros(paddingL,1);
for i = 1:signalTTLNum
    spacerTTL(round(1+(fs*signalTTLiti*(i-1))):round(fs*(ttlDur+signalTTLiti*(i-1)))) = 1;
end

zeroSig = zeros(paddingL,1);
signalVector = [zeroSig,spacerTTL];

%recalculate ramp times for onset and offset in samples
onRampDur = onRampDur*fs; 
offRampDur = offRampDur*fs;
remainingPoints = L-onRampDur-offRampDur;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
%generate the overall ramp profile
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;

%% Calculate laser lag and incorporate into duration of sound file

%first, compensate for loss of time due to second TTL pulse
laserLag = laserLag - 0.004;

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

%% Generating sound files

%calculate amplitude of sounds
toneAmplA = 10^-((maxDB-toneDBA)/20);
toneAmplB = 10^-((maxDB-toneDBB)/20);

%generate actual sounds for A
toneWaveA = sin(2*pi*(toneFreqA/fs)*(1:paddingL))';
toneWaveA = toneWaveA.*rampProfile;
toneWaveA = toneWaveA * toneAmplA;
paddedWaveA = zeros(paddingL,1);
paddedWaveA(1:size(toneWaveA,1)) = toneWaveA;
soundVectorA = [paddedWaveA,ttlSig];
soundVectorALaser = [paddedWaveA,laserTTL];

%generate actual sounds for B
toneWaveB = sin(2*pi*(toneFreqB/fs)*(1:paddingL))';
toneWaveB = toneWaveB.*rampProfile;
toneWaveB = toneWaveB * toneAmplB;
paddedWaveB = zeros(paddingL,1);
paddedWaveB(1:size(toneWaveB,1)) = toneWaveB;
soundVectorB = [paddedWaveB,ttlSig];

%% Now execute delivery of sounds.

laserSwitch = 0; %switch that controls what kind of trial to play

for i = 1:numTrials
    %check if you are at a boundary. If so, then execute the margin
    %signal
    %Also executes a pause command
    if i == 2*baselineReps + 1 | i == 2*baselineReps + 2*pairingReps + 1
        pause(master(i,1)/2)
        sound(signalVector,fs)
        pause(master(i,1)/2)
        if i == 2*baselineReps + 1
            laserSwitch = 1;
            disp('Ending Baseline Period, Starting Laser Period')
        elseif i == 2*baselineReps + 2*pairingReps + 1
            laserSwitch = 0;
            disp('Ending Laser Period, Starting Recovery Period')
        end
    else
        pause(master(i,1))
    end
    %determine which tone to play
    if master(i,2) == 1
        if laserSwitch == 1
            sound(soundVectorALaser,fs)
            disp(strcat('Trial',num2str(i),'/',num2str(numTrials),'Freq:',num2str(master(i,3)),'DB:',num2str(master(i,4))))
        elseif laserSwitch == 0
            sound(soundVectorA,fs)
            disp(strcat('Trial',num2str(i),'/',num2str(numTrials),'Freq:',num2str(master(i,3)),'DB:',num2str(master(i,4))))
        end
    elseif master(i,2) == 0
        sound(soundVectorB,fs)
        disp(strcat('Trial',num2str(i),'/',num2str(numTrials),'Freq:',num2str(master(i,3)),'DB:',num2str(master(i,4))))
    end
    
    
end

%save the information! (important!!)
soundData = struct;
soundData.Delays = master(:,1);
soundData.BaselineRepetitions = baselineReps;
soundData.LaserRepetitions = pairingReps;
soundData.FinishRepetitions = finishReps;
soundData.Frequencies = master(:,3);
soundData.dBs = master(:,4);
soundData.Master = master;
soundData.ToneDuration = toneDur;
soundData.TargetFreq = toneFreqA;
soundData.TargetDB = toneDBA;
soundData.FoilFreq = toneFreqB;
soundData.FoilDB = toneDBB;
soundData.LaserLag = laserLag;


save(fullfile(pname,fname),'soundData');






