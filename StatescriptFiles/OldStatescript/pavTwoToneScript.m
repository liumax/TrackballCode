%This is the Matlab Script

global scQtUserData;
          
%% UI prompt:
prompt = {'Mouse ID:',...  
    'Weight:',...  
    'Big Reward (msec):',...
    'Small Reward (msec):',...       
    'Trials (even #):',...          
    'Tone Duration (msec):',...
    'Tone Amplitude (dB, max 100):',...
    'Reward Delay (msec):',...
    'Big Tone (Hz):',...       
    'Small Tone (Hz):',...       
    'ITI (msec):',...
    'ITI (Longest) (msec):',...
    'sessionID:',...        
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'LickTask:';
num_lines=1;
def={'','','200','50','200','500','100','3000','','','10000','50000','1',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

%% store all the shits
i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.bigRew = str2num(answer{i});i=i+1;
scQtUserData.smallRew = str2num(answer{i});i=i+1;
scQtUserData.totalTrials = str2num(answer{i});i=i+1;
scQtUserData.soundDur = str2num(answer{i});i=i+1;
scQtUserData.soundAmp = str2num(answer{i});i=i+1;
scQtUserData.rewDelay = str2num(answer{i});i=i+1;
scQtUserData.bigTone = str2num(answer{i});i=i+1;
scQtUserData.smallTone = str2num(answer{i});i=i+1;
scQtUserData.ITI = str2num(answer{i});i=i+1;
scQtUserData.ITIRange = str2num(answer{i});i=i+1;
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;
scQtUserData.taskID = 'twoTonePavlovComputer';

%% now lets start calculations. 
%calculate ITIs, use exponential distribution
k = 2.5;
p = (1-exp(-k))*rand(scQtUserData.totalTrials,1);
tau = (scQtUserData.ITIRange-scQtUserData.ITI)/k;
x = round(scQtUserData.ITI-6000 + (-log(1-p))*tau); 
scQtUserData.Master(:,1) = x;

%calculate lag for reward delivery
k = 2.5;
p = (1-exp(-k))*rand(scQtUserData.totalTrials,1);
tau = (300)/k; %adjusted 170808 to try and improve behavior. 
x = round(scQtUserData.rewDelay + (-log(1-p))*tau); 
scQtUserData.RewDelayMatrix = x;

%determine when to present tones!
randInd = randperm(scQtUserData.totalTrials);
randInd(randInd<=scQtUserData.totalTrials/2) = 1; %1 will represent small rewards
randInd(randInd>scQtUserData.totalTrials/2) = 2; %2 will represent large rewards
scQtUserData.Master(:,2) = randInd;
%determine rewSize order
scQtUserData.Master(:,3) = zeros(scQtUserData.totalTrials,1);
scQtUserData.Master(randInd == 1,3) = scQtUserData.smallRew;
scQtUserData.Master(randInd == 2,3) = scQtUserData.bigRew;
%the last thing is to have a lick window such that you enforce a no lick
%period before the cue delivery. 
scQtUserData.lickWindow = 2000; %generates a 2 second window for licking. This will only be useful in the full behavior and not in training.

%now I need to do all the prep for the sounds

%FIXED PARAMETERS
fs = 192000; %sampling frequency in Hz

calibChart = [4000	2	0.7943282347
4287.09385	0	1
4594.79342	0.5	0.9440608763
4924.577653	1	0.8912509381
5278.031643	1.5	0.8413951416
5656.854249	1.5	0.8413951416
6062.866266	2.2	0.7762471166
6498.019171	1.5	0.8413951416
6964.404506	2	0.7943282347
7464.263932	4.5	0.5956621435
8000	6.2	0.4897788194
8574.1877	4.3	0.6095368972
9189.58684	3.6	0.660693448
9849.155307	6.3	0.4841723676
10556.06329	5.1	0.5559042573
11313.7085	3.8	0.645654229
12125.73253	1.5	0.8413951416
12996.03834	3.5	0.6683439176
13928.80901	3.7	0.6531305526
14928.52786	3.8	0.645654229
16000	3.5	0.6683439176
17148.3754	2.5	0.7498942093
18379.17368	2	0.7943282347
19698.31061	6.2	0.4897788194
21112.12657	7.8	0.4073802778
22627.417	8.75	0.3651741273
24251.46506	10	0.316227766
25992.07668	14	0.1995262315
27857.61803	12	0.2511886432
29857.05573	13.2	0.2187761624
32000	15.6	0.1659586907
34296.7508	16.5	0.1496235656
36758.34736	18.2	0.1230268771
39396.62123	20	0.1
42224.25314	16	0.1584893192
45254.834	18.7	0.1161448614
48502.93013	14.2	0.19498446
51984.15337	15.5	0.1678804018
55715.23605	11.7	0.2600159563
59714.11146	11.7	0.2600159563
64000	10	0.316227766];

bigToneAmp = interp1(calibChart(:,1),calibChart(:,3),scQtUserData.bigTone);
smallToneAmp = interp1(calibChart(:,1),calibChart(:,3),scQtUserData.smallTone);

if isnan(bigToneAmp)
    error('bigToneAmp ISNAN')
elseif isnan(smallToneAmp)
    error('smallToneAmp ISNAN')
end

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
waveSmall = (sin(2*pi*(scQtUserData.smallTone/fs)*(1:L))') .* rampProfile;

%calculate amplitude
toneDB = 10^-((100-scQtUserData.soundAmp)/20);

paddedWave = zeros(paddingL,1);
paddedWave(1:size(waveBig,1)) = waveBig;
soundBig = [paddedWave*toneDB*bigToneAmp,ttlSig];
paddedWave = zeros(paddingL,1);
paddedWave(1:size(waveSmall,1)) = waveSmall;
soundSmall = [paddedWave*toneDB*smallToneAmp,ttlSig];

scQtUserData.ToneBig = soundBig;
scQtUserData.ToneSmall = soundSmall;


%store information about time/date
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.tripSwitch = 0;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
scQtUserData.failTrig = 0; %trigger for failure of the graphing code. 

pause(0.5);

%% Display these to statescript log file for storage. 
sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''weight:', scQtUserData.weight,''')']);
sendScQtControlMessage(['disp(''bigReward:', num2str(scQtUserData.bigRew),''')']);
sendScQtControlMessage(['disp(''smallReward:', num2str(scQtUserData.smallRew),''')']);
sendScQtControlMessage(['disp(''totalTrials:', num2str(scQtUserData.totalTrials),''')']);
sendScQtControlMessage(['disp(''soundDur:', num2str(scQtUserData.soundDur),''')']);
sendScQtControlMessage(['disp(''soundAmp:', num2str(scQtUserData.soundAmp),''')']);
sendScQtControlMessage(['disp(''rewDelay:', num2str(scQtUserData.rewDelay),''')']);
sendScQtControlMessage(['disp(''bigTone:', num2str(scQtUserData.bigTone),''')']);
sendScQtControlMessage(['disp(''smallTone:', num2str(scQtUserData.smallTone),''')']);
sendScQtControlMessage(['disp(''ITIShort:', num2str(scQtUserData.ITI),''')']);
sendScQtControlMessage(['disp(''ITILong:', num2str(scQtUserData.ITIRange),''')']);
sendScQtControlMessage(['disp(''lickWindow:', num2str(scQtUserData.lickWindow),''')']);
sendScQtControlMessage(['disp(''taskID:', scQtUserData.taskID,''')']);
sendScQtControlMessage(['disp(''date:', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''time:', scQtUserData.time,''')']);
sendScQtControlMessage(['disp(''sessionID:', scQtUserData.sessionID,''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

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
scQtUserData.licks = zeros(1000,4);
scQtUserData.lickCounter = 1;

scQtUserData.lickHist = zeros(80,2); %This is optimized for looking at an 8 second window
%with 2 sec before sound onset, sound, and 3 seconds after. Set for 100 ms
%bins. First column is for small reward, second column for big rewards.

scQtUserData.lickAxes = [-2:0.1:5.9]; %axis for histogram
%send initial information to the mbed
sendScQtControlMessage(['lickWind =',num2str(scQtUserData.lickWindow)]);
sendScQtControlMessage(['toneRewDel =',num2str(scQtUserData.RewDelayMatrix(1))]);
sendScQtControlMessage(['trPhase = 0']);
sendScQtControlMessage(['signalDel =3000']); %this is the delay after reward delivery before triggering next thing. 
sendScQtControlMessage(['disp(''StartSession'')']);
