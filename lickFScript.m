%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  1
    'Min Reward (msec):',...       2
    'Max Reward (msec):',...       
    'Blocks:',...           5
    'Block Size:',...       6
    'Sound Duration (ms)',...    7
    'Weight:',...           8
    'sessionID:',...        9
    'Notes:'}; %the bracket is to end the prompt     10
dlg_title = 'LickTask:';
num_lines=1;
def={'','50','200','8','50','3000','','1',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

% Connect to Microcontroller and store MetaData:
% setComVal; % Call script to set com port string for this computer
% sHandle = scConnect(comValStr,@pokeE_Callback);
% pause(1);

waterWindow=1000;

i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.minRew = str2num(answer{i});i=i+1;
scQtUserData.maxRew = str2num(answer{i});i=i+1;
scQtUserData.blocks = str2num(answer{i});i=i+1;
scQtUserData.blockSize = str2num(answer{i});i=i+1;
scQtUserData.soundDur = str2num(answer{i});i=i+1;
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.taskID = 'LickTask';
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;

scQtUserData.minITI=scQtUserData.soundDur+6000;
scQtUserData.maxITI=scQtUserData.minITI+5000;

scQtUserData.timeDelay = 4000; %generates 3 second delay so callback triggered later.
scQtUserData.lickWindow = 2000; %generates a 2 second window for licking. this will be used before sound
%so that there is an enforced no lick period, and after sound onset, for
%lick-triggering of reward.

scQtUserData.waterWindow=waterWindow;
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.tripSwitch = 0;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
% scQtUserData.trInit = false;
% scQtUserData.trBait = false;

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''minRew:', num2str(scQtUserData.minRew),''')']);
sendScQtControlMessage(['disp(''maxRew:', num2str(scQtUserData.maxRew),''')']);
sendScQtControlMessage(['disp(''blocks:', num2str(scQtUserData.blocks),''')']);
sendScQtControlMessage(['disp(''blockSize:', num2str(scQtUserData.blockSize),''')']);
sendScQtControlMessage(['disp(''minITI:', num2str(scQtUserData.minITI),''')']);
sendScQtControlMessage(['disp(''maxITI:', num2str(scQtUserData.maxITI),''')']);
sendScQtControlMessage(['disp(''waterWindow:', num2str(scQtUserData.waterWindow),''')']);
sendScQtControlMessage(['disp(''soundDur:', num2str(scQtUserData.soundDur),''')']);
sendScQtControlMessage(['disp(''weight:', scQtUserData.weight,''')']);
sendScQtControlMessage(['disp(''taskID:', scQtUserData.taskID,''')']);
sendScQtControlMessage(['disp(''date:', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''time:', scQtUserData.time,''')']);
sendScQtControlMessage(['disp(''sessionID:', scQtUserData.sessionID,''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

pause(1) %Need to put all my timings in before this stuff

triallength=scQtUserData.blocks*scQtUserData.blockSize;

%master array for all calculations
master=zeros(triallength,11);

%master(:,1) determines reward size. fills in rewards sizes so can pull by
%trial number
masterRand = rand(1);
if masterRand>0.5
    master(:,1)=scQtUserData.minRew;
    for x=2:2:scQtUserData.blocks
        master(1+(x-1)*scQtUserData.blockSize:x*scQtUserData.blockSize,1)=scQtUserData.maxRew;
    end
else
    master(:,1)=scQtUserData.maxRew;
    for x=2:2:scQtUserData.blocks
        master(1+(x-1)*scQtUserData.blockSize:x*scQtUserData.blockSize,1)=scQtUserData.minRew;
    end
end

%master(:,2) will implement an exponential for the ITI distribution, with
%random noise inserted.

k = 2.5;
p = (1-exp(-k))*rand(triallength,1);
tau = (scQtUserData.maxITI-scQtUserData.minITI)/k;
% x = round(scQtUserData.minITI - scQtUserData.soundDur +
% (-log(1-p))*tau)-2000; THIS MAY BE PROBLEM WITH ITIs
x = round(scQtUserData.minITI + (-log(1-p))*tau)-scQtUserData.lickWindow-scQtUserData.timeDelay; 
%This is adjusted to allow for SoundOff to be trigger for matlab. -1000 is
%adjustment for pre delays (1000 each) in callback. Post delays do not affect 
%when 'SoundOff' appears, so they are not accounted for in here. 
%%6/18/2015 edit, changing this to 4000 to account for predelays and
%%timedelay, which is necessary to keep scripts functional. 
master(:,2)= x;

%master(:,3) will calculate delay time from tone presentation. Will use
%flat distribution utilizing random numbers.

delayRatio=waterWindow;
master(:,3)=round(scQtUserData.soundDur-(rand(triallength,1)*delayRatio));

%This next line subtracts master(:,3) from master(:,2). This is because the
%statescript nests these times, so I must substract them properly to get
%proper ITIs.

master(:,2)=master(:,2)-master(:,3);

%master(:,4) calculates if reward is delivered; 1 means delivery, 0 means
%none. 6/18/2015 deleted so rewards delivered every time. 

%master(:,5) determines probability of laser; 1 means delivery, 0 means
%none. this has been removed, since there is no plan currently for lasers

%This will be for correctly completed trials
master(:,6)=zeros(triallength,1);

%This is for failed trials with no licks
master(:,7)=zeros(triallength,1);

%This is for trials with badlicks
master(:,8)=zeros(triallength,1);

%This is for Sound Times (triggered by SoundOn)
master(:,10)=zeros(triallength,1);

%This is for calculation of ITIs 
master(:,11)=zeros(triallength,1);

scQtUserData.master=master;

scQtUserData.licks = zeros(1000,1);
scQtUserData.lickCounter = 1;

scQtUserData.lickHist = zeros(80,2); %This is optimized for looking at an 8 second window
%with 2 sec before sound onset, sound, and 3 seconds after. Set for 100 ms
%bins. First column is for small reward, second column for big rewards.

scQtUserData.lickAxes = [-2:0.1:5.9];

sendScQtControlMessage(['soundDur=',num2str(scQtUserData.soundDur)]);
sendScQtControlMessage(['lickWindow=',num2str(scQtUserData.lickWindow)]);
sendScQtControlMessage(['timeDelay=',num2str(scQtUserData.timeDelay)]);
sendScQtControlMessage(['disp(''StartSession'')']);
