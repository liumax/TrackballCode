%This is the Matlab Script

%I HAVE A HARD OFFSET FOR TIMING. This is due to a roughly 200-350ms delay
%between execution of sound(s,sf) and actual delivery of sound, as detected
%by TTL pulse. 
hardOffset = 350; %hard offset in ms.

global scQtUserData;
          
% UI prompt:
prompt = {'File Name:',...  
    'Tone Duration (s):',...
    'Tone Repetitions:',...       
    'Tone ITI (min, ms):',...        
    'Tone ITI (max, ms):',...
    'Starting Frequency (Hz):',...
    'Ending Frequency (Hz):',...
    'Frequency Jumps (octaves):',...
    'Starting DB (dB):',...         
    'Ending DB (dB):',...    
    'DB Jumps (dBs):',...
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'basicTuningCurve:';
num_lines=1;
def={'','0.1','20','500','1000','4000','64000','1','100','60','20',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

% Connect to Microcontroller and store MetaData:
% setComVal; % Call script to set com port string for this computer
% sHandle = scConnect(comValStr,@pokeE_Callback);
% pause(1);

i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.fileName = answer{i};i=i+1;
toneDur = str2num(answer{i});i=i+1;
toneReps = str2num(answer{i});i=i+1;
postPauseMin = str2num(answer{i});i=i+1;
postPauseMax = str2num(answer{i});i=i+1;
startF = str2num(answer{i});i=i+1;
endF = str2num(answer{i});i=i+1;
octFrac = str2num(answer{i});i=i+1;
startdB = str2num(answer{i});i=i+1;
enddB = str2num(answer{i});i=i+1;
dBSteps = str2num(answer{i});i=i+1;
scQtUserData.notes = answer{i};i=i+1;
scQtUserData.taskID = 'basicTuningCurve';

% safety check: confirm neither frequency is lower than 4000 Hz.
if startF < 4000 | endF < 4000
    error('FREQUENCY TOO LOW. TRIAL ABORTED')
end

%hardcoded values
prePause = 100; %pre-pause in ms.
fs = 192000;
L = toneDur * fs; %duration of signal.
paddingL = L + fs*0.2; %duration with safety padding
maxdB = 100; %maximum dB output. 
onRampDur = 0.005*fs;
offRampDur = 0.005*fs;
%store some into global
scQtUserData.paddingL = paddingL;
scQtUserData.L = L;
scQtUserData.fs = 192000;
scQtUserData.ToneDur = toneDur*1000;

%Date/Time
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
scQtUserData.tripSwitch = 0; %thing to shut things down

pause(0.2);

%% calculate all the things I'll need for the tuning curve:

octRange = log2(endF/startF);
freqs = zeros(octRange/octFrac+1,1);
freqs(1) = startF;
for i = 2:octRange/octFrac+1
    freqs (i) = freqs(i-1)*(2^octFrac);
end

%this generates a vector of all decibel steps
dBs = [startdB:-dBSteps:enddB];
amps = ones(length(dBs),1);
for i = 1:length(amps)
    amps(i) = amps(i)*10^-((maxdB-dBs(i))/20);
end

%list of all frequency/dB pairs
fullList = zeros(length(freqs)*length(dBs),3);
counter = 1;

%generates a large array with all combinations of frequencies and dBs/amps
for i = 1:length(freqs)
    for j = 1:length(dBs)
        fullList(counter:counter+toneReps-1,1) = freqs(i);
        fullList(counter:counter+toneReps-1,2) = dBs(j);
        fullList(counter:counter+toneReps-1,3) = amps(j);
        counter = counter+toneReps;
    end
end

% generates a randomized set of indices for calling things from fullList.
fillIndex = zeros(toneReps*length(freqs)*length(dBs),1);
counter = 1;
for i = 1:toneReps*length(dBs)
    fillIndex(counter:counter+length(freqs)-1) = randperm(length(freqs));
    counter = counter + length(freqs);
end

%creates array for storing everything. Also fills index with information
%from above.
master = zeros(length(freqs)*length(dBs)*toneReps,3);

for i = 1:length(fillIndex)
    master(i,1) = freqs(fillIndex(i)); %fills in frequency based on fill index
    holder = find(fullList(:,1) == master(i,1)); %finds where frequency matches with fullIndex
    holderSize = size(holder,1); %finds how many match, generates random integer within that range.
    [x randIndex] = sort(rand(holderSize,1));
    master(i,2) = fullList(holder(randIndex(1)),2); %fills in dB and amplitude
    master(i,3) = fullList(holder(randIndex(1)),3); %amplitude value
    fullList(holder(randIndex(1)),:) = []; %removes from fullList so that you sample without replacement.
end

%save total number of trials, store in global
totalTrials = length(master);
scQtUserData.TotalTrials = totalTrials;

%calculates ITIs
k = 2.5;
p = (1-exp(-k))*rand(totalTrials,1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 
itiTime = x;

master(:,4) = x+hardOffset;

%save master to global
scQtUserData.Master = master;

%ramp times for onset and offset in seconds
remainingPoints = L-onRampDur-offRampDur;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;
scQtUserData.rampProfile = rampProfile;

%this makes the profile for the TTL signal
ttlSig = zeros(paddingL,1);
ttlSig(1:2*fs/1000) = 1;
scQtUserData.ttlSig = ttlSig;



%% save file with information.
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


save(strcat(scQtUserData.fileName,'.mat'),'soundData');

%% 
sendScQtControlMessage(['disp(''Filename:', scQtUserData.fileName,''')']);
sendScQtControlMessage(['disp(''ToneDuration:', num2str(toneDur),''')']);
sendScQtControlMessage(['disp(''ToneRepetitions:', num2str(toneReps),''')']);
sendScQtControlMessage(['disp(''MinITI:', num2str(postPauseMin),''')']);
sendScQtControlMessage(['disp(''MaxITI:', num2str(postPauseMax),''')']);
sendScQtControlMessage(['disp(''StartFrequency:', num2str(startF),''')']);
sendScQtControlMessage(['disp(''EndFrequency:', num2str(endF),''')']);
sendScQtControlMessage(['disp(''OctaveJump:', num2str(octFrac),''')']);
sendScQtControlMessage(['disp(''StartDB:', num2str(startdB),''')']);
sendScQtControlMessage(['disp(''EndDB:', num2str(enddB),''')']);
sendScQtControlMessage(['disp(''dbSteps:', num2str(dBSteps),''')']);
sendScQtControlMessage(['disp(''TotalTrials:', num2str(totalTrials),''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

pause(1) %Need to put all my timings in before this stuff

sendScQtControlMessage(['disp(''StartSession'')']);
sendScQtControlMessage(['toneDur = ',num2str(scQtUserData.ToneDur)]); 
