%overall structure of this system is to basically provide pseudorandomized
%delivery of tones, with variable delays between laser onset and tone
%onset. This is meant to allow me to determine which timings are most
%effective at altering response properties. 

%open dialog box for saving file
[fname pname] = uiputfile('test1.mat');


%% User Controlled Parameters

%tone duration
toneDur = 0.1; %tone duration in seconds
onRampDur = 0.005; %both onramp and offramp will be cosine waves
offRampDur = 0.005;

%tone amplitude and frequency
toneFreq = 16000; %if want white noise, enter 0
toneDB = 80;

%task structure (assumes laser is delivered over 400ms)
repsPerLag = 100; %number of tone presentations per block
laserLags = [0 0.5 0.75 1 1.25]; %different positive lags in seconds (tone follows laser). 
%NOTE: 0 lag is no laser. 
postPauseMin = 5; %max and minimum pause times. In seconds.
postPauseMax = 7;

%% Fixed Parameters
fs = 192000; %sampling frequency in Hz
prePause = 0.1; %pause in seconds before tone delivery
maxDB = 100;
ttlDur = 0.002;
ttlITI = 0.004; %iti between TTL onsets for laser.

%% Calculations!

%calculate amplitude
toneAmp = 10^-((maxDB-toneDB)/20);

%determine length of sound file
L = round((max(laserLags)+toneDur)*fs); %want to generate file with maximum file size to cover all laser possibilities
paddingL = round(L*1.1); %make the padding equal to 10% of overall sound file

%determine ITI
k = 2.5;
p = (1-exp(-k))*rand(length(laserLags)*repsPerLag,1);
tau = (postPauseMax-postPauseMin)/k;
x = (postPauseMin + (-log(1-p))*tau); 

%determine length of set of frequencies to be played
listLength = length(laserLags);

%calculate pseudorandomized order
listDesig = zeros(listLength*repsPerLag,1);
counter = 1;
for i = 1:repsPerLag
    listDesig(counter:counter+listLength-1) = randperm(listLength);
    counter = counter + listLength;
end

%generate master array for storage of information
master = zeros(length(listDesig),5);
master(:,1) = toneFreq;
master(:,2) = toneDB;
master(:,3) = toneAmp;
master(:,4) = x;
master(:,5) = listDesig;

%calculate tone wave
if toneFreq == 0
    y = wgn(L,1,0);
    %filters white gaussian noise (courtesy of RYAN MORRILL)
    Wn = 4e3/(0.5*fs); % pass above 3.9 kHz 160121 Adjusted to 4kHz
    n = 1000; % 1000th order filter (slower? but 100-order was too low)
    b = fir1(n, Wn, 'high'); 
    toneWave = filtfilt(b,1,y);
elseif toneFreq >=4000
    toneWave = sin(2*pi*(toneFreq/fs)*(1:L))';
elseif toneFreq < 4000 & toneFreq ~= 0
    error('Tone Frequency Too Low')
end

%recalculate ramp times for onset and offset in samples
onRampDur = onRampDur*fs; 
offRampDur = offRampDur*fs;
remainingPoints = L-onRampDur-offRampDur;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
rampProfile(1:max(laserLags)*fs) = 0;
rampProfile(max(laserLags)*fs:max(laserLags)*fs+onRampDur) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;

%generate sound vector (since I dont ever need to recalculate!)
finalWave = toneWave.*rampProfile*toneAmp;
paddedWave = zeros(paddingL,1);
paddedWave(1:size(finalWave,1)) = finalWave;

%check to confirm ramp falls within expected limits
if length(find(rampProfile==1)) > fs*toneDur;
    error('Failure to generate correct profile, Excessively long')
end

%generate a set of TTLs to represent different laser lags.
%first, make a generic TTL for tone. 
ttlSig = zeros(paddingL,1);
ttlSig(max(laserLags)*fs:max(laserLags)*fs+ttlDur*fs) = 1;
%next, make set of TTLs for different lags
ttlHolder = zeros(paddingL,listLength);
for i = 1:listLength
    if laserLags(i) == 0
        ttlHolder(:,i) = ttlSig;
        ttlHolder(1:round(1+ttlDur*fs),i) = 1;
    elseif laserLags(i) == max(laserLags)
        ttlHolder(:,i) = ttlSig;
        ttlHolder(1:round(1+ttlDur*fs),i) = 1;
        ttlHolder(round(ttlITI*fs):round((ttlITI+ttlDur)*fs),i) = 1;
    else
        %find difference (how far back to move)
        diffInd = round((max(laserLags)-laserLags(i))*fs);
        ttlHolder(:,i) = ttlSig;
        ttlHolder(diffInd:round(diffInd+ttlDur*fs),i) = 1;
        ttlHolder(round(diffInd+ttlITI*fs):round(diffInd+(ttlITI+ttlDur)*fs),i) = 1;
    end
end

%% Perform Sound Presentations

for i = 1:length(master)
    pause(prePause)
    soundVector = [paddedWave,ttlHolder(:,master(i,5))];
    sound(soundVector,fs);
    disp(strcat('Trial:',num2str(i),'/',num2str(length(master)),'Lag:',num2str(laserLags(master(i,5)))))
    pause(master(i,4));
end

soundData = struct;
soundData.Master = master;
soundData.ToneDuration = toneDur;
soundData.ToneRepetitions = repsPerLag;
soundData.Frequencies = master(:,1);
soundData.dBs = master(:,2);
soundData.Amplitudes = master(:,3);
soundData.Delays = master(:,4);
soundData.LaserDesig = master(:,5);
soundData.LaserLags = laserLags;

save(fullfile(pname,fname),'soundData');

