%this is code to generate a white with amplitude modulation

fs = 192000;       %sampling frequency
tsec = 0.1;             %length of sound in seconds
L = tsec*fs;        % length of sound in samples
A = 1;                  % Amplitude of the envelope
freq = 21000; %frequency in Hz of tone.
ramp = 0.005; %ramp duration (on and off, in seconds)

%for unmodulated
envelope = ones(1,L);

%generates linear ramps of specified duration
ramps = ones(L,1);
ramps(1:(ramp*fs)) = [0:1/(ramp*fs):1-1/(ramp*fs)];
ramps(end-(ramp*fs):end) = [1:-1/(ramp*fs):0];

%combines envelopes
envelope = ramps .* envelope';

%generates Tone!
toneWave = sin(2*pi*(freq/fs)*(1:L))';

%combines things together
soundWave = toneWave.*envelope;
finalSound = zeros(length(soundWave),2);
finalSound(:,1) = soundWave;
finalSound(:,2) = soundWave;

audioplayer(soundWave,fs,16,5)

fileName = '10kHzTest.wav'
audiowrite(fileName,finalSound,fs);

