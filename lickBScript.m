%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  1
    'Min Reward (msec):',...       2
    'Max Reward (msec):',...       3
    'Reward Prob:',...      4
    'Reverse Contingency?:',...       5
    'Blocks:',...           6
    'Block Size:',...       7
    'Sound:'...             8
    'Sound Duration (ms)',...    9
    'Weight:',...           10
    'sessionID:',...        11
    'Notes:'}; %the bracket is to end the prompt     12
dlg_title = 'LickTask:';
num_lines=1;
def={'','50','200','1','n','8','50','','3000','','1',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

% Connect to Microcontroller and store MetaData:
% setComVal; % Call script to set com port string for this computer
% sHandle = scConnect(comValStr,@pokeE_Callback);
% pause(1);

waterWindow=1000;
baitDur=1000;

i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.minRew = str2double(answer{i});i=i+1;
scQtUserData.maxRew = str2double(answer{i});i=i+1;
scQtUserData.rewProb = str2double(answer{i});i=i+1;
scQtUserData.reversal = answer{i};i=i+1;
scQtUserData.blocks = str2double(answer{i});i=i+1;
scQtUserData.blockSize = str2double(answer{i});i=i+1;
scQtUserData.sound = answer{i};i=i+1;
scQtUserData.soundDur = str2double(answer{i});i=i+1;
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.taskID = 'LickTask';
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;

scQtUserData.minITI=scQtUserData.soundDur+4000;
scQtUserData.maxITI=scQtUserData.minITI+5000;

scQtUserData.waterWindow=waterWindow;
scQtUserData.baitDur=baitDur;
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.tripSwitch = 0;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
% scQtUserData.trInit = false;
% scQtUserData.trBait = false;

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID: ', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''minRew: ', num2str(scQtUserData.minRew),''')']);
sendScQtControlMessage(['disp(''maxRew: ', num2str(scQtUserData.maxRew),''')']);
sendScQtControlMessage(['disp(''rewProb: ', num2str(scQtUserData.rewProb),''')']);
sendScQtControlMessage(['disp(''reversal: ', scQtUserData.reversal,''')']);
sendScQtControlMessage(['disp(''blocks: ', num2str(scQtUserData.blocks),''')']);
sendScQtControlMessage(['disp(''blockSize: ', num2str(scQtUserData.blockSize),''')']);
sendScQtControlMessage(['disp(''minITI: ', num2str(scQtUserData.minITI),''')']);
sendScQtControlMessage(['disp(''maxITI: ', num2str(scQtUserData.maxITI),''')']);
sendScQtControlMessage(['disp(''waterWindow: ', num2str(scQtUserData.waterWindow),''')']);
sendScQtControlMessage(['disp(''sound: ', scQtUserData.sound,''')']);
sendScQtControlMessage(['disp(''soundDur: ', num2str(scQtUserData.soundDur),''')']);
sendScQtControlMessage(['disp(''weight: ', scQtUserData.weight,''')']);
sendScQtControlMessage(['disp(''taskID: ', scQtUserData.taskID,''')']);
sendScQtControlMessage(['disp(''date: ', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''time: ', scQtUserData.time,''')']);
sendScQtControlMessage(['disp(''sessionID: ', scQtUserData.sessionID,''')']);
sendScQtControlMessage(['disp(''notes: ', scQtUserData.notes,''')']);
sendScQtControlMessage(['disp(''baitDur: ', num2str(scQtUserData.baitDur),''')']);

pause(1) %Need to put all my timings in before this stuff

triallength=scQtUserData.blocks*scQtUserData.blockSize;

%master array for all calculations
master=zeros(triallength,12);

%master(:,1) determines reward size. fills in rewards sizes so can pull by
%trial number

if ~isempty(strfind(scQtUserData.sound,'lo'))
    master(:,1)=scQtUserData.minRew;
elseif ~isempty(strfind(scQtUserData.sound,'hi'))
    master(:,1)=scQtUserData.maxRew;
elseif ~isempty(strfind(scQtUserData.sound,'both'))
    starter=rand;
    if starter > 0.5
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
end

%master(:,2) will implement an exponential for the ITI distribution, with
%random noise inserted.

k = 2.5;
p = (1-exp(-k))*rand(triallength,1);
tau = (scQtUserData.maxITI-scQtUserData.minITI)/k;
% x = round(scQtUserData.minITI - scQtUserData.soundDur +
% (-log(1-p))*tau)-2000; THIS MAY BE PROBLEM WITH ITIs
x = round(scQtUserData.minITI + (-log(1-p))*tau)-1000; 
%This is adjusted to allow for SoundOff to be trigger for matlab. -1000 is
%adjustment for pre delays (1000 each) in callback. Post delays do not affect 
%when 'SoundOff' appears, so they are not accounted for in here. 
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
%none.

master(:,4)=rand(triallength,1);
master(master(:,4)>=scQtUserData.rewProb,4)=0;
master(master(:,4)<scQtUserData.rewProb,4)=1;

%master(:,5) determines probability of laser; 1 means delivery, 0 means
%none.

% % master(:,5)=rand(triallength,1);
% % master(master(:,5)>=scQtUserData.rewProb,5)=0; PROBABLY BUGGY
% % master(master(:,5)<scQtUserData.rewProb,5)=1;

%This will be for pre-cue licks
master(:,6)=zeros(triallength,1);

%This is for anticipatory licks
master(:,7)=zeros(triallength,1);

%This is for consummatory licks
master(:,8)=zeros(triallength,1);

%This is for licks in all other intervals
master(:,9)=zeros(triallength,1);

%This is found Sound Times (triggered by SoundOff)
master(:,10)=zeros(triallength,1);

%This is for calculation of ITIs (marked by sound-off)
master(:,11)=zeros(triallength,1);

%This is for calculation of which sound to deliver
if ~isempty(strfind(scQtUserData.reversal,'y'))
    if ~isempty(strfind(scQtUserData.sound,'lo'))
        master(:,12) = 2;
    elseif ~isempty(strfind(scQtUserData.sound,'hi'))
        master(:,12) = 1;
    elseif ~isempty(strfind(scQtUserData.sound,'both'))
        master(:,12) = master(:,1);
        bigRew=max(master(:,12));
        smallRew=min(master(:,12));
        master(master(:,12)==bigRew,12)=1;
        master(master(:,12)==smallRew,12)=2;
    else
        disp 'NO SOUND SELECTED ABORT'
        pause
    end
elseif ~isempty(strfind(scQtUserData.reversal,'n'))
    if ~isempty(strfind(scQtUserData.sound,'lo'))
        master(:,12) = 1;
    elseif ~isempty(strfind(scQtUserData.sound,'hi'))
        master(:,12) = 2;
    elseif ~isempty(strfind(scQtUserData.sound,'both'))
        master(:,12) = master(:,1);
        bigRew=max(master(:,12));
        smallRew=min(master(:,12));
        master(master(:,12)==bigRew,12)=2;
        master(master(:,12)==smallRew,12)=1;
    else
        disp 'NO SOUND SELECTED ABORT'
        pause
    end
end

scQtUserData.master=master;

sendScQtControlMessage(['soundDur=',num2str(scQtUserData.soundDur)]);
sendScQtControlMessage(['baitDur=',num2str(baitDur)]);
sendScQtControlMessage(['disp(''SoundOff'')']);