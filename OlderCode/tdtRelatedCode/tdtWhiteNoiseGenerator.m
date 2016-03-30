%control over settings here
toneReps = 100; %number of repetitions of each tone/amplitude pair
toneDur = 0.1; %tone duration in seconds
ttlDur = 0.01; %duration of signaling TTL in seconds
fs = 192000; %sampling frequency in Hz
L = toneDur*fs;

prePause = 0.2; %pause in seconds before tone
postPause = 0.4; %pause in seconds after tone

warningCheck = (postPause - toneDur)<0;
if warningCheck == 1
    disp('TONE DURATION LONGER THAN ITI')
end

%generate white gaussian noise of correct size
y = wgn(fs,1,0);

%code to plot fft. validated this does appropriate filtering on 2016-01-15
% testy = fft(y);
% P2 = abs(testy/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% plot(f,P1)

%filters white gaussian noise (courtesy of RYAN MORRILL)
Wn = 3.9e3/(0.5*fs); % pass above 3.9 kHz
n = 1000; % 1000th order filter (slower? but 100-order was too low)
b = fir1(n, Wn, 'high'); 
audio_data = filtfilt(b,1,y);
    
%code to plot fft. validated this does appropriate filtering on 2016-01-15
% hold on
% tester = fft(audio_data);
% P2 = abs(tester/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% plot(f,P1,'r')

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
    disp(tonereps - i)
    pause(postPause)
end
