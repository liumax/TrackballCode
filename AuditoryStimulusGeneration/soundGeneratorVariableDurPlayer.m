%this code is to generate a target tone of set intensity at a variety of
%tone durations. 

[fname pname] = uiputfile('durTester.mat');

%% Parameters
%Tone Parameters
toneFreq = 4000; %tone frequency in Hz
% toneFreq = 'white'; %option for white noise!
toneDB = 100; %tone loudness in dB

%Presentation Parameters
toneReps = 100; %repetitions of each tone length
toneDur = [0.1 1 5]; %tone duration in seconds
minITI = max(toneDur) * 1.5; %minimum ITI between tone onsets
maxITI = max(toneDur) * 3; %maximum ITI between tone onsets

%fixed parameters
ttlDur = 0.002; %duration of TTL pulse marking tone
fs = 192000; %sampling frequency of sound card
maxDB = 100; %maximum output.
prePause = 0.1; %prepause in seconds for loop so that we dont get crashes
timeBuffer = 0.5; %time buffer in seconds from end of tone. This is to avoid tone cutoffs.
rampDur = 0.005; %ramp duration in seconds. Will be cosine ramp. 

%% Safety Checks
if toneFreq < 4000;
    error('Insufficient Tone Frequency, will damage speaker')
end

%% Calculations
%calculate amplitude
toneAmp= 10^-((maxDB-toneDB)/20);

%calculate number of total repetitions
totalReps = toneReps * length(toneDur);

%generates random ITIs using exponential function
k = 2.5;
p = (1-exp(-k))*rand(totalReps,1);
tau = (maxITI-minITI)/k;
pauseTimes = minITI + (-log(1-p))*tau; 

%calculate longest possible signal
L = round(fs*max(toneDur));
paddingL = L + round(fs*timeBuffer);

%generate TTL signals
genericSig = zeros(paddingL,1);
genericSig(1:round(ttlDur*fs)) = 1;

ttlSig = cell(length(toneDur),1);

for i = 1:length(toneDur)
    ttlSig{i} = genericSig;
    ttlSig{i}(round(toneDur(i)*fs):round((toneDur(i)+ttlDur)*fs)) = 1;
end

%generate ramp profiles
sampleRampDur = rampDur*fs;
onRampProfile = (cos((0:1:sampleRampDur)/sampleRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:sampleRampDur)/sampleRampDur*pi)+1)/2;

genericProf = zeros(paddingL,1);
genericProf(1:length(onRampProfile)) = onRampProfile;

rampProfs = cell(length(toneDur),1);

for i = 1:length(toneDur)
    rampProfs{i} = genericProf;
    rampProfs{i}(length(onRampProfile):round(toneDur(i)*fs)) = 1;
    rampProfs{i}(round(toneDur(i)*fs)-length(offRampProfile)+1:round(toneDur(i)*fs)) = offRampProfile;
end

%generate sine wave/white noise
if isnumeric(toneFreq)
    toneWave=sin(2*pi*(toneFreq/fs)*(1:paddingL))';
elseif strfind(toneFreq,'white')
    y = wgn(paddingL,1,0);
    %filters white gaussian noise (courtesy of RYAN MORRILL)
    Wn = 4e3/(0.5*fs); % pass above 3.9 kHz 160121 Adjusted to 4kHz
    n = 1000; % 1000th order filter (slower? but 100-order was too low)
    %160121 Adjusted from 1000 to 100000, this gets greater attenuation, which is
    %necessary for my higher volumes. This achieves 90dB attenuation for all
    %frequencies below 4kHz.
    %160121 Returned to 1000. Increasing n to 100000 causes failure of
    %filtfilt. Checked by doing fft(X,n) with n = number of samples, seems to
    %filter out everything below 4kHz very effectively. 
    b = fir1(n, Wn, 'high'); 
    toneWave = filtfilt(b,1,y);
end


%combine these things!
finalVectors = cell(length(toneDur),1);
for i = 1:length(toneDur)
    finalVectors{i} = [(toneWave*toneAmp).*rampProfs{i},ttlSig{i}];
end

%generate order of tones
toneOrder = zeros(totalReps,1);
toneCounter = 1;
for i =1:toneReps
    toneOrder(toneCounter:toneCounter+length(toneDur)-1) = randperm(length(toneDur));
    toneCounter = toneCounter + length(toneDur);
end



%% Execute Code

for i = 1:totalReps
    pause(prePause)
    sound(finalVectors{toneOrder(i)},fs)
    pause(pauseTimes(i))
end



soundData = struct;
soundData.ToneFreq = toneFreq;
soundData.ToneDB = toneDB;
soundData.ToneAmplitude = toneAmp;
soundData.ToneDurations =toneDur;
soundData.TTLDur = ttlDur;
soundData.RampDur = rampDur;

save(fullfile(pname,fname),'soundData');

