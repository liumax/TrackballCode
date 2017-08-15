
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
octFrac = 0.2; %fractions of octaves to move
%White noise?
whiteNoise = 1; %1 for white noise inclusion.

%dB
startdB = 100; %starting decibel value
enddB = 60; %lowest decibel value
dbSteps = 10; %resolution of decible steps

%confirm pauses are longer than double the tone length
warningCheck = (postPauseMin/1000 - 2*toneDur)<0;
if warningCheck == 1
    error('TONE DURATION LONGER THAN ITI')
end

%% Speaker Calibration (external #1)

calibChart = [4000	2	0.7943282347
4287.09385	0	1
4594.79342	0.5	0.9440608763
4924.577653	1	0.8912509381
5278.031643	1.5	0.8413951416
5656.854249	1.5	0.8413951416
6062.866266	2.2	0.7762471166
6498.019171	1.5	0.8413951416
6964.404506	2	0.7943282347
7464.263932	4.5	0.5956621435
8000	6.2	0.4897788194
8574.1877	4.3	0.6095368972
9189.58684	3.6	0.660693448
9849.155307	6.3	0.4841723676
10556.06329	5.1	0.5559042573
11313.7085	3.8	0.645654229
12125.73253	1.5	0.8413951416
12996.03834	3.5	0.6683439176
13928.80901	3.7	0.6531305526
14928.52786	3.8	0.645654229
16000	3.5	0.6683439176
17148.3754	2.5	0.7498942093
18379.17368	2	0.7943282347
19698.31061	6.2	0.4897788194
21112.12657	7.8	0.4073802778
22627.417	8.75	0.3651741273
24251.46506	10	0.316227766
25992.07668	14	0.1995262315
27857.61803	12	0.2511886432
29857.05573	13.2	0.2187761624
32000	15.6	0.1659586907
34296.7508	16.5	0.1496235656
36758.34736	18.2	0.1230268771
39396.62123	20	0.1
42224.25314	16	0.1584893192
45254.834	18.7	0.1161448614
48502.93013	14.2	0.19498446
51984.15337	15.5	0.1678804018
55715.23605	11.7	0.2600159563
59714.11146	11.7	0.2600159563
64000	10	0.316227766];

%% FIXED PARAMETERS
fs = 192000; %sampling frequency in Hz
maxdB = 100; %maximum decibel output

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
if whiteNoise == 1
    freqs(end+1) = 0;
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

%170815 generate interpolation for db and amplitude adjustments
dbChange = interp1(calibChart(:,1),calibChart(:,2),freqs);
ampChange = interp1(calibChart(:,1),calibChart(:,3),freqs);

for i = 1:length(freqs)
    for j = 1:length(dBs)
        fullList(counter,1) = freqs(i);
        %find the frequency
        if freqs(i) >0
            
            fullList(counter,2) = dBs(j) - dbChange(i);
            fullList(counter,3) = amps(j) * ampChange(i);
        else
            fullList(counter,2) = dBs(j);
            fullList(counter,3) = amps(j);
        end
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

%generate white noise to be used repeatedly.
%generate white gaussian noise of correct size
y = wgn(fs*toneDur,1,0);

%filters white gaussian noise (courtesy of RYAN MORRILL)
Wn = 4e3/(0.5*fs); % pass above 3.9 kHz 160121 Adjusted to 4kHz
n = 1000; % 1000th order filter (slower? but 100-order was too low)
b = fir1(n, Wn, 'high'); 
audio_data = filtfilt(b,1,y);


%% Perform actual Tuning
for i = 1:length(master)
    pause(prePause)
    toneFreq = master(i,1);
    if toneFreq == 0;
        toneWave = audio_data;
    else
        toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
    end
    toneAmpl = master(i,3);
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
soundData.WhiteNoise = whiteNoise;



save(fullfile(pname,fname),'soundData');
