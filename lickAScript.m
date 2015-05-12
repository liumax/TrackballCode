%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  1
    'Min Reward (msec):',...       2
    'Max Reward (msec):',...       3
    'Reward Prob:',...      4
    'Laser Prob:',...       5
    'Blocks:',...           6
    'Block Size:',...       7
    'Sound Duration (ms)',...    8
    'Training? 1 = yes',... 
    'Weight:',...           9
    'sessionID:',...        10
    'Notes:'}; %the bracket is to end the prompt     11
dlg_title = 'LickTask:';
num_lines=1;
def={'','50','200','1','0','8','50','3000','','','1',''};
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
scQtUserData.laserProb = str2double(answer{i});i=i+1;
scQtUserData.blocks = str2double(answer{i});i=i+1;
scQtUserData.blockSize = str2double(answer{i});i=i+1;
scQtUserData.soundDur = str2double(answer{i});i=i+1;
scQtUserData.training = str2double(answer{i});i=i+1;
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


% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
% scQtUserData.trInit = false;
% scQtUserData.trBait = false;

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID: ', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''minRew: ', num2str(scQtUserData.minRew),''')']);
sendScQtControlMessage(['disp(''maxRew: ', num2str(scQtUserData.maxRew),''')']);
sendScQtControlMessage(['disp(''rewProb: ', num2str(scQtUserData.rewProb),''')']);
sendScQtControlMessage(['disp(''laserProb: ', num2str(scQtUserData.rewProb),''')']);
sendScQtControlMessage(['disp(''blocks: ', num2str(scQtUserData.blocks),''')']);
sendScQtControlMessage(['disp(''blockSize: ', num2str(scQtUserData.blockSize),''')']);
sendScQtControlMessage(['disp(''minITI: ', num2str(scQtUserData.minITI),''')']);
sendScQtControlMessage(['disp(''maxITI: ', num2str(scQtUserData.maxITI),''')']);
sendScQtControlMessage(['disp(''waterWindow: ', num2str(scQtUserData.waterWindow),''')']);
sendScQtControlMessage(['disp(''soundDur: ', num2str(scQtUserData.soundDur),''')']);
sendScQtControlMessage(['disp(''Training: ', num2str(scQtUserData.training),''')']);
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
master=zeros(triallength,10);

%master(:,1) determines reward size. fills in rewards sizes so can pull by
%trial number
master(:,1)=scQtUserData.minRew;
for x=2:2:scQtUserData.blocks
    master(1+(x-1)*scQtUserData.blockSize:x*scQtUserData.blockSize,1)=scQtUserData.maxRew;
end

%master(:,2) will implement an exponential for the ITI distribution, with
%random noise inserted.

k = 2.5;
p = (1-exp(-k))*rand(triallength,1);
tau = (scQtUserData.maxITI-scQtUserData.minITI)/k;
x = round(scQtUserData.minITI - scQtUserData.soundDur + (-log(1-p))*tau); %This is adjusted to allow for SoundOff to be trigger for matlab
master(:,2)= x;

%master(:,3) will calculate delay time from tone presentation. Will use
%flat distribution utilizing random numbers.

delayRatio=waterWindow;
master(:,3)=round(scQtUserData.soundDur-(rand(triallength,1)*delayRatio));

%master(:,4) calculates if reward is delivered; 1 means delivery, 0 means
%none.

master(:,4)=rand(triallength,1);
master(master(:,4)>=scQtUserData.rewProb,5)=0;
master(master(:,4)<scQtUserData.rewProb,5)=1;

%master(:,5) determines probability of laser; 1 means delivery, 0 means
%none.

master(:,5)=rand(triallength,1);
master(master(:,5)>=scQtUserData.rewProb,5)=0;
master(master(:,5)<scQtUserData.rewProb,5)=1;

%This will be for pre-cue licks
master(:,6)=zeros(triallength,1);

%This is for anticipatory licks
master(:,7)=zeros(triallength,1);

%This is for consummatory licks
master(:,8)=zeros(triallength,1);

scQtUserData.master=master;

sendScQtControlMessage(['soundDur=',num2str(scQtUserData.soundDur)]);
sendScQtControlMessage(['baitDur=',num2str(baitDur)]);
sendScQtControlMessage('SoundOff');
