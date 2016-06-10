%this is code to generate a white with amplitude modulation

fs = 44000;       %sampling frequency
tsec = 0.1;             %length of sound in seconds
L = tsec*fs;        % length of sound in samples
A = 1;                  % Amplitude of the envelope
m = 0.5               % modulation depth of the envelope
Fm = 100;     %frequency of modulation in Hz
ramp = 0.005; %ramp duration (on and off, in seconds)

%for unmodulated
envelope = ones(1,L);

%for modulated
% envelope = A * [1 + m *sin(2*pi*Fm/fs*(1:L) - pi/2)]; %equation for the envelope

%generates linear ramps of specified duration
ramps = ones(L,1);
ramps(1:(ramp*fs)) = [0:1/(ramp*fs):1-1/(ramp*fs)];
ramps(end-(ramp*fs):end) = [1:-1/(ramp*fs):0];

%combines envelopes
envelope = ramps .* envelope';

%generates white noise
sigma = 2;      % standard deviation of gaussian noise
mu = 0;           % mean of white gaussian noise

noise = sigma*randn(L,1)+mu;
noise = noise/max(abs(noise));

%combines things together
soundWave = noise.*envelope;
finalSound = zeros(length(soundWave),2);
finalSound(:,1) = soundWave;
finalSound(:,2) = soundWave;

fileName = '002_noise100ms.wav'
audiowrite(fileName,finalSound,fs);

