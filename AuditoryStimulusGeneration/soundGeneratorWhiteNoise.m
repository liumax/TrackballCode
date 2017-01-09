%this code is to generate white noise. have not added laser functionality.

%control over settings here
toneReps = 100; %number of repetitions of each tone/amplitude pair
totalReps = toneReps;
toneDur = 0.1; %tone duration in seconds
ttlDur = 0.01; %duration of signaling TTL in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs;

%pausing times!
minPause = 1;
maxPause = 2;

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

%ramp times for onset and offset in seconds
onRampDur = 0.005; 
offRampDur = 0.005;

%this code generates linear ramps for onset and offset. this reduces issues
%with sharp white noises
rampProfile = ones(L,1);
rampProfile(1:(onRampDur*fs)) = [0:1/(onRampDur*fs):1-1/(onRampDur*fs)];
rampProfile(end-(onRampDur*fs):end) = [1:-1/(onRampDur*fs):0];

%this makes the profile for the TTL signal
ttlSig = zeros(L,1);
ttlSig(1:fs*ttlDur) = 1;

%adds ramps to white noise
finalWave = audio_data.*rampProfile;

%makes two column matrix for sound and TTL output
soundVector = [finalWave,ttlSig];

for i = 1:toneReps
    pause(prePause)
    sound(soundVector,fs);
    disp(toneReps - i)
    pause(x(i))
end
