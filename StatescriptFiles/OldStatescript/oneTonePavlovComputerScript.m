%This is the Matlab Script

global scQtUserData;
          
%% UI prompt:
prompt = {'Mouse ID:',...  
    'Weight:',...  
    'Reward (msec):',...    
    'Trials (even #):',...          
    'Tone Duration (msec):',...
    'Tone Amplitude (dB, max 100):',...
    'Reward Delay (msec):',...
    'Tone (Hz):',...            
    'ITI (msec):',...
    'ITI (Longest) (msec):',...
    'sessionID:',...        
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'LickTask:';
num_lines=1;
def={'','','200','200','500','100','3000','','10000','50000','1',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

%% store all the shits
i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.bigRew = str2num(answer{i});i=i+1;
scQtUserData.totalTrials = str2num(answer{i});i=i+1;
scQtUserData.soundDur = str2num(answer{i});i=i+1;
scQtUserData.soundAmp = str2num(answer{i});i=i+1;
scQtUserData.rewDelay = str2num(answer{i});i=i+1;
scQtUserData.bigTone = str2num(answer{i});i=i+1;
scQtUserData.ITI = str2num(answer{i});i=i+1;
scQtUserData.ITIRange = str2num(answer{i});i=i+1;
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;
scQtUserData.taskID = 'oneTonePavlovComputer';

%% now lets start calculations. 
%calculate ITIs, use exponential distribution
k = 2.5;
p = (1-exp(-k))*rand(scQtUserData.totalTrials,1);
tau = (scQtUserData.ITIRange-scQtUserData.ITI)/k;
x = round(scQtUserData.ITI-6000 + (-log(1-p))*tau); 

scQtUserData.Master(:,1) = x;
scQtUserData.Master(:,2) = ones(scQtUserData.totalTrials,1);
%determine rewSize order
scQtUserData.Master(:,3) = ones(scQtUserData.totalTrials,1)*scQtUserData.bigRew;
%the last thing is to have a lick window such that you enforce a no lick
%period before the cue delivery. 
scQtUserData.lickWindow = 2000; %generates a 2 second window for licking. This will only be useful in the full behavior and not in training.

%now I need to do all the prep for the sounds

%FIXED PARAMETERS
fs = 192000; %sampling frequency in Hz

%calculations for tone
%determine length of file based on sound card sampling rate
L = scQtUserData.soundDur/1000*fs; %number of samples at correct sampling frequency
paddingL = round(L*1.5); %adds 50% time as buffer

%recalculate ramp times for onset and offset in samples
onRampDur = 0.005*fs; 
offRampDur = 0.005*fs;
remainingPoints = L-onRampDur-offRampDur;
onRampProfile = (cos((0:1:onRampDur)/onRampDur*pi-pi)+1)/2;
offRampProfile = (cos((0:1:offRampDur)/offRampDur*pi)+1)/2;
rampProfile = ones(L,1);
rampProfile(1:onRampDur+1) = onRampProfile;
rampProfile(end-offRampDur:end) = offRampProfile;
%this makes the profile for the TTL signal
ttlSig = zeros(paddingL,1);
ttlSig(1:5*fs/1000) = 1;
%generate the two sounds
waveBig = (sin(2*pi*(scQtUserData.bigTone/fs)*(1:L))') .* rampProfile;
% waveSmall = (sin(2*pi*(scQtUserData.smallTone/fs)*(1:L))') .* rampProfile;

%calculate amplitude
toneDB = 10^-((100-scQtUserData.soundAmp)/20);

paddedWave = zeros(paddingL,1);
paddedWave(1:size(waveBig,1)) = waveBig;
soundBig = [paddedWave*toneDB,ttlSig];
% paddedWave = zeros(paddingL,1);
% paddedWave(1:size(waveSmall,1)) = waveSmall;
% soundSmall = [paddedWave*toneDB,ttlSig];

scQtUserData.ToneBig = soundBig;
% scQtUserData.ToneSmall = soundSmall;


%store information about time/date
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.tripSwitch = 0;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number

pause(0.5);

%% Display these to statescript log file for storage. 
sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''weight:', scQtUserData.weight,''')']);
sendScQtControlMessage(['disp(''bigReward:', num2str(scQtUserData.bigRew),''')']);
% sendScQtControlMessage(['disp(''smallReward:', num2str(scQtUserData.smallRew),''')']);
sendScQtControlMessage(['disp(''totalTrials:', num2str(scQtUserData.totalTrials),''')']);
sendScQtControlMessage(['disp(''soundDur:', num2str(scQtUserData.soundDur),''')']);
sendScQtControlMessage(['disp(''soundAmp:', num2str(scQtUserData.soundAmp),''')']);
sendScQtControlMessage(['disp(''rewDelay:', num2str(scQtUserData.rewDelay),''')']);
sendScQtControlMessage(['disp(''bigTone:', num2str(scQtUserData.bigTone),''')']);
% sendScQtControlMessage(['disp(''smallTone:', num2str(scQtUserData.smallTone),''')']);
sendScQtControlMessage(['disp(''ITIShort:', num2str(scQtUserData.ITI),''')']);
sendScQtControlMessage(['disp(''ITILong:', num2str(scQtUserData.ITIRange),''')']);
sendScQtControlMessage(['disp(''lickWindow:', num2str(scQtUserData.lickWindow),''')']);
sendScQtControlMessage(['disp(''taskID:', scQtUserData.taskID,''')']);
sendScQtControlMessage(['disp(''date:', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''time:', scQtUserData.time,''')']);
sendScQtControlMessage(['disp(''sessionID:', scQtUserData.sessionID,''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);
sendScQtControlMessage(['disp(''trPhase = 0'')']);
pause(1) %Need to put all my timings in before this stuff

%% generate space in structure for storage of information that I care about!
%This is for soundOn times. This allows for calculations of all the other
%things!
scQtUserData.soundOn = zeros(scQtUserData.totalTrials,1);
%this is for licking latency
scQtUserData.lickLat = zeros(scQtUserData.totalTrials,1);


%variables for tracking licking. Each entry here will be a combination of
%the lick time relative to the sound, the trial number, and the type of
%trial.
scQtUserData.licks = zeros(100000,3);
scQtUserData.lickCounter = 1;
scQtUserData.cueTime = zeros(scQtUserData.totalTrials,1);
scQtUserData.firstLick = zeros(scQtUserData.totalTrials,1);

scQtUserData.lickHist = zeros(80,2); %This is optimized for looking at an 8 second window
%with 2 sec before sound onset, sound, and 3 seconds after. Set for 100 ms
%bins. First column is for small reward, second column for big rewards.

scQtUserData.lickAxes = [-2:0.1:5.9]; %axis for histogram
%send initial information to the mbed
sendScQtControlMessage(['lickWind =',num2str(scQtUserData.lickWindow)]);
sendScQtControlMessage(['toneRewDel =',num2str(scQtUserData.rewDelay)]);
sendScQtControlMessage(['trPhase = 0']);
sendScQtControlMessage(['signalDel =3000']); %this is the delay after reward delivery before triggering next thing. 
sendScQtControlMessage(['disp(''StartSession'')']);
