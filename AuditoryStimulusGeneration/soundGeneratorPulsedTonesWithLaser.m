%this code is to generate pulsed sounds, either white noise or tones.
%Sounds will be played as repeated short tones in blocks, and these blocks
%will be spaced by a larger ITI

[fname pname] = uiputfile('test1.mat');

%%control over settings here
%general settings
fs = 192000; %sampling frequency in Hz
prePause = 0.1; %pause in seconds before tone delivery

%control over timing of blocks of tone pips
totalReps = 10; %number of repetitions of tone blocks
totalITI = 10; %ITI of tone blocks in seconds. DOES NOT INCLUDE TONE DUR

%control over parameters of individual tone pips
% toneFreq = 8000; %either frequency in Hz or 'white'
toneFreq = 'white';
toneDB = 60;%amplitude in dB with 100 set as max with driver at no attenuation, speakers set to max
maxDB = 100;%DONT ADJUST
toneAmp= 10^-((maxDB-toneDB)/20);

%control for tone pips 
controlFreq = 8000;
% controlFreq = 'white';
controlDB = 100;
controlAmp = 10^-((maxDB-controlDB)/20);

indivReps = 15; %number of tone pips per presentation
indivITI = 0.2; %ITI of tone pips in seconds. This includes the tone duration.
toneDur = 0.1; %tone pip duration in seconds
onRampDur = 0.005; %duration in seconds of on-ramp. This will be a cosine ramp
offRampDur = 0.005;%duration in seconds of off-ramp. This will be a cosine ramp

blockDur = (indivReps)*indivITI*fs*1.1; %time of block in samples, additional 10% fudge factor to avoid cutoff.
tonePipDur = toneDur*fs; %duration of tone pip in samples. Will be useful in generating template for sound to insert into larger tone file.
indivITIdur = indivITI*fs; %iti in samples

toneTTL = 0.002; %duration in seconds of ttl that marks the tone onset

%laser controls
laserSwitch = 1; %1 for laser, 0 for no laser
laserTiming = 5; %pick the tone pip that the laser is going to be paired with
laserTTL = 0.002; %duration of laser TTL
laserTTLiti = 0.006; %ITI for laser TTLs. recall this is from start of TTL to start of TTL
laserTTLnum = 2; %number of laser TTLs

%%safety checks
if toneDur > indivITI
    error('Tone Duration Shorter than ITI')
end

if onRampDur + offRampDur > toneDur
    error('Tone is too short for ramps')
end

if ~ischar(toneFreq) & toneFreq < 4000;
    error('Tone Frequency is Unsafely Low')
end

if laserTTLiti - laserTTL < 0.002
    error('Laser ITI is insufficiently spaced')
end

if maxDB - toneDB < 0
    error('DB value greater than maximum')
end

if laserSwitch ~= 1 & laserSwitch ~= 0
    error('Laser Mode is Invalid')
end

if laserTiming > indivReps
    error('Laser Timing is Invalid')
end

%%start prepping sound pip file. This code makes the small tone pip sound and TTL.

%first, make the amplitude envelope
rampProfile = ones(tonePipDur,1);
rampProfile(1:(onRampDur*fs)) = [0:1/(onRampDur*fs):1-1/(onRampDur*fs)];
rampProfile(end-(onRampDur*fs):end) = [1:-1/(onRampDur*fs):0];
rampProfile = rampProfile * toneAmp; %adjusts for amplitude

%this makes the profile for the TTL signal
pipTTL = zeros(tonePipDur,1);
pipTTL(1:fs*toneTTL) = 1;

laserPipTTL = zeros(tonePipDur,1);
laserPipTTL(1:fs*laserTTL) = 1;
laserPipTTL(laserTTLiti*fs:fs*(laserTTL+laserTTLiti)) = 1;

%fill tone space with white noise or tone
if ischar(toneFreq) & strfind(toneFreq,'white')
    %generate white gaussian noise of correct size
    y = wgn(fs*toneDur,1,0);
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
    audio_data = filtfilt(b,1,y);
    tonePip = audio_data.*rampProfile;
elseif ~ischar(toneFreq)
    tonePip = sin(2*pi*(toneFreq/fs)*(1:tonePipDur))';
    tonePip = tonePip .* rampProfile;
end

%fill control tone space with white noise or tone
if ischar(controlFreq) & strfind(controlFreq,'white')
    %generate white gaussian noise of correct size
    y = wgn(fs*toneDur,1,0);
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
    audio_data = filtfilt(b,1,y);
    controlPip = audio_data.*rampProfile;
elseif ~ischar(controlFreq)
    controlPip = sin(2*pi*(controlFreq/fs)*(1:tonePipDur))';
    controlPip = controlPip .* rampProfile;
end

%%now put the small pips together! This will also incorporate laser. This
%%generates a two column matrix that incorporates sound (column 1) and ttl
%%(column 2)

soundBlock = zeros(blockDur,2);
startTimes = 1:indivITIdur:indivITIdur*indivReps;

%fills in the array!
for i = 1:indivReps
    soundBlock(startTimes(i):startTimes(i)-1+tonePipDur,1) = tonePip;
    if laserSwitch & i == laserTiming
        soundBlock(startTimes(i):startTimes(i)-1+tonePipDur,2) = laserPipTTL;
    else
        soundBlock(startTimes(i):startTimes(i)-1+tonePipDur,2) = pipTTL;
    end
end

%generates a laser free sound array
controlBlock = zeros(blockDur,2);

for i = 1:indivReps
    controlBlock(startTimes(i):startTimes(i)-1+tonePipDur,1) = controlPip;
    controlBlock(startTimes(i):startTimes(i)-1+tonePipDur,2) = pipTTL;
end

%generate actual pause time (this is because sound(s,sf) is instantaneously
%executed, and then moves onto the next line of code.
pauseDur = totalITI + (indivReps)*indivITI;

%tone execution
if laserSwitch == 1; 
    for i = 1:totalReps*2
        pause(prePause)
        if mod(i,2) == 0 %finds evens, plays laser
            soundVector = soundBlock;
            sound(soundVector,fs);
        else %plays control on odds
            soundVector = controlBlock;
            sound(soundVector,fs);
        end
        totalReps*2 - i
        totalReps*2
        pause(pauseDur);
    end
else
    for i = 1:totalReps
        pause(prePause)
        soundVector = controlBlock;
        sound(soundVector,fs);
        totalReps - i
        totalReps
        pause(pauseDur);
    end
end

soundData = struct;
soundData.TargetFrequency = toneFreq;
soundData.TargetDB = toneDB;
soundData.TargetAmplitude = toneAmp;
soundData.ControlFrequency = controlFreq;
soundData.ControlDB = controlDB;
soundData.ControlAmplitude = controlAmp;

soundData.BlockRepetitions = totalReps;
soundData.BlockITI = totalITI;
soundData.PipRepetitions = indivReps;
soundData.PipITI = indivITI;
soundData.OnRampDur = onRampDur;
soundData.OffRampDur = offRampDur;
soundData.LaserOn = laserSwitch;
soundData.LaserTiming = laserTiming;
soundData.LaserTTL = laserTTL;
soundData.LaserTTLITI = laserTTLiti;

save(fullfile(pname,fname),'soundData');



