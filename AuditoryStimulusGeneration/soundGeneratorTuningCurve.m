
%opens dialogue box to change saved file name
[fname pname] = uiputfile('test1.mat');

%% USER ADJUSTED PARAMETERS
%Duration and Repetition
toneReps = 30; %number of repetitions of each tone/amplitude pair
toneDur = 0.1; %tone duration in seconds
onRampDur = 0.005;
offRampDur = 0.005;

%Timing
prePause = 0.1; %pause in seconds before tone
postPauseMin = 500; %pause in milliseconds after tone
postPauseMax = 1000; %pause in milliseconds after tone

%Frequency
startF = 4000; %starting frequency in Hz
endF = 64000; %ending frequency in Hz
octFrac = 0.5; %fractions of octaves to move

%dB
maxdB = 100; %maximum decibel output
startdB = 100; %starting decibel value
enddB = 60; %lowest decibel value
dbSteps = 20; %resolution of decible steps

%confirm pauses are longer than double the tone length
warningCheck = (postPauseMin/1000 - 2*toneDur)<0;
if warningCheck == 1
    error('TONE DURATION LONGER THAN ITI')
end

%% FIXED PARAMETERS
fs = 192000; %sampling frequency in Hz

%% Calculations for Tuning Curve
%determine length of file based on sound card sampling rate
L = toneDur*fs; %number of samples at correct sampling frequency
paddingL = L + fs*0.1; %adds 0.1 seconds of padding to the end of the tone to ensure things are not cut off.

%Calculate the range of frequencies that will be used. 
octRange = log2(endF/startF);
freqs = zeros(octRange/octFrac+1,1);
freqs(1) = startF;
for i = 2:octRange/octFrac+1
    freqs (i) = freqs(i-1)*(2^octFrac);
end

%Calculate range of all dBs that will be used.
dBs = [startdB:-dbSteps:enddB];
amps = ones(length(dBs),1);
for i = 1:length(amps)
    amps(i) = amps(i)*10^-((maxdB-dBs(i))/20);
end

%List the full set of combinations of frequency and dB. This should produce
%all combinations, BUT NOT all tone repetitions
fullList = zeros(length(freqs)*length(dBs),3);
counter = 1;

for i = 1:length(freqs)
    for j = 1:length(dBs)
        fullList(counter,1) = freqs(i);
        fullList(counter,2) = dBs(j);
        fullList(counter,3) = amps(j);
        counter = counter+1;
    end
end
%Compute length of the list to do a check against hard calculation
listLength = length(fullList);
%check this is correct
if listLength ~= length(freqs)*length(dBs)
    error('ListLength Incorrect')
end

%generate a n x 1 vector that indicates all indices for every tuning trial.
%This generates a pseudorandom set of indices for determining the frequency
%and dB of the target trial.
listDesig = zeros(listLength*toneReps,1);
counter = 1;
for i = 1:toneReps
    listDesig(counter:counter + listLength - 1) = randperm(listLength);
    counter = counter + listLength;
end

%creates array for storing everything. Also fills index with information
%from above.
master = zeros(length(freqs)*length(dBs)*toneReps,4);
master(:,1:3) = fullList(listDesig,1:3);

%generates pauses with exponential function
k = 2.5;
p = (1-exp(-k))*rand(length(master),1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 

master(:,4) = x/1000;

%recalculate ramp times for onset and offset in samples
onRampDur = onRampDur*fs; 
offRampDur = offRampDur*fs;
remainingPoints = L-onRampDur-offRampDur;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;


%this makes the profile for the TTL signal
ttlSig = zeros(paddingL,1);
ttlSig(1:2*fs/1000) = 1;


%% Perform actual Tuning
for i = 1:length(master)
    pause(prePause)
    toneFreq = master(i,1);
    toneAmpl = master(i,3);
    toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
    finalWave = toneWave.*rampProfile;
    finalWave = finalWave*toneAmpl;
    paddedWave = zeros(paddingL,1);
    paddedWave(1:size(finalWave,1)) = finalWave;
    soundVector = [paddedWave,ttlSig];
    sound(soundVector,fs);
    disp(strcat('Trial:',num2str(i),'/',num2str(length(master)),' Frequency:',num2str(master(i,1)),' DB:',num2str(master(i,2))))
    pause(master(i,4))
end


%generate matrix with trial data.

trialMatrix = zeros(length(master),4);
trialMatrix(:,1) = [1:1:length(master)];
trialMatrix(:,2:4) = master(:,1:3);

trialHeaders = ['Trial','Frequency','dB','Amplitude'];

soundData = struct;
soundData.Delays = master(:,4);
soundData.ToneRepetitions = toneReps;
soundData.ToneDuration = toneDur;
soundData.TrialMatrix = trialMatrix;
% soundData.UniqueFreqs = freqs;
% soundData.UniqueDBs = dBs;
soundData.Frequencies = master(:,1);
soundData.dBs = master(:,2);
soundData.Amplitudes = master(:,3);



save(fullfile(pname,fname),'soundData');
