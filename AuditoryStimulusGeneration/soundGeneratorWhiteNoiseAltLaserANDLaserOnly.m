%this code is to generate white noise with alternating laser
[fname pname] = uiputfile('test1.mat');
%control over settings here
toneReps = 100; %number of repetitions of each condition
totalReps = toneReps*3;
toneDur = 0.2; %tone duration in seconds
ttlDur = 0.002; %duration of signaling TTL in seconds
fs = 192000; %sampling frequency in Hz
L = (toneDur*2)*fs;

%time before tone that I want the laser to come on in seconds
laserLag = 0.4;
laserJitter = 0.1;
%find the number of samples to pad.

laserITI = 0.006; %time betwen laser pulses

jitterVec = (rand(toneReps,1)-0.5)*laserJitter+laserLag;

%pausing times!
minPause = 3;
maxPause = 5;

prePause = 0.1; %pause in seconds before tone

%generates random ITIs using exponential function
k = 2.5;
p = (1-exp(-k))*rand(totalReps,1);
tau = (maxPause-minPause)/k;
x = minPause + (-log(1-p))*tau; 


warningCheck = (minPause - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

%generate pseudorandom order
orderList = zeros(totalReps,1);
orderHold = 1;
for i = 1:toneReps
    orderList(orderHold:orderHold + 2) = randperm(3);
    orderHold = orderHold + 3;
end

%generate white gaussian noise of correct size
y = wgn(L,1,0);

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

%ramp times for onset and offset in seconds
onRampDur = 0.005; 
offRampDur = 0.005;

%this code generates linear ramps for onset and offset. this reduces issues
%with sharp white noises
rampProfile = ones(L,1);
rampProfile(1:(onRampDur*fs)) = [0:1/(onRampDur*fs):1-1/(onRampDur*fs)];
rampProfile(round(toneDur*fs):round((toneDur+offRampDur)*fs)) = [1:-1/(onRampDur*fs):0];
rampProfile(round((toneDur+offRampDur)*fs):end) = 0;

%this makes the profile for the TTL signal for tone
ttlSig = zeros(L,1);
ttlSig(1:fs*ttlDur) = 1;

%make ttl sig for laser only. 
laserOnlyTTL = zeros(L,1);
laserOnlyTTL(1:fs*ttlDur) = 1;
laserOnlyTTL(round(laserITI*fs):round(laserITI*fs+fs*ttlDur)) = 1;

laserOnlyVector = [zeros(L,1),laserOnlyTTL];

%adds ramps to white noise
finalWave = audio_data.*rampProfile;


%makes two column matrix for sound and TTL output
soundVector = [finalWave,ttlSig];


for i = 1:toneReps*2
    pause(prePause)
    if orderList(i) == 1 %play just the sound
        sound(soundVector,fs)
    elseif orderList(i) == 2 %play sound and laser
        laserLadd = round(jitterVec(i/2) * fs);
        laserTTLSig = zeros(L+laserLadd,1);
        laserTTLSig(laserLadd+1:end) = ttlSig;
        laserTTLSig(1:round(fs*ttlDur)) = 1;
        laserTTLSig(round(laserITI*fs):round(laserITI*fs+fs*ttlDur)) = 1;
        finalLaserWave = audio_data.*rampProfile;
        finalLaserWave(laserLadd+1:end+laserLadd) = finalLaserWave;
        finalLaserWave(1:laserLadd) = 0;
        soundLaserVector = [finalLaserWave,laserTTLSig];
        sound(soundLaserVector,fs)
    elseif orderList(i) == 3 %play just laser
        sound(laserOnlyVector,fs)
    end
    disp(totalReps - i)
    pause(x(i))
end

soundData = struct;
soundData.ToneDuration = toneDur;
soundData.ToneRepetitions = toneReps;
soundData.ToneDesignation = orderList;
soundData.OrderTags = {'NoiseOnly','Noise+Laser','LaserOnly'};
soundData.LaserLag = laserLag;
soundData.ITI = x;

save(fullfile(pname,fname),'soundData');