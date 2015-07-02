%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  1
    'Reward (msec):',...       2
    'Trials (even #):',...           5
    'Reward Tone (p/s):',...       6
    'Weight:',...           8
    'sessionID:',...        9
    'Notes:'}; %the bracket is to end the prompt     10
dlg_title = 'LickTask:';
num_lines=1;
def={'','50','200','','','1',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

% Connect to Microcontroller and store MetaData:
% setComVal; % Call script to set com port string for this computer
% sHandle = scConnect(comValStr,@pokeE_Callback);
% pause(1);

i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.rewSize = str2num(answer{i});i=i+1;
scQtUserData.totalTrials = str2num(answer{i});i=i+1;
scQtUserData.rewTone = answer{i};i=i+1;
scQtUserData.soundDur = 3000;
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

scQtUserData.waterWindow=1000; %window starting from end of sound where reward can be delivered.
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.tripSwitch = 0;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
% scQtUserData.trInit = false;
% scQtUserData.trBait = false;

if (~isempty(strfind(scQtUserData.rewTone,'p')))
    scQtUserData.toneRule = 1;
elseif (~isempty(strfind(scQtUserData.rewTone,'s')))
    scQtUserData.toneRule = 0;
else
    disp('FAILED PARSING OF SOUND RULE')
    pause
end

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''rewSize:', num2str(scQtUserData.rewSize),''')']);
sendScQtControlMessage(['disp(''totalTrials:', num2str(scQtUserData.totalTrials),''')']);
sendScQtControlMessage(['disp(''rewTone:', num2str(scQtUserData.rewTone),''')']);
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

%calculate the sound to deliver
random1 = randperm(scQtUserData.totalTrials);
random1(random1<=scQtUserData.totalTrials/2) = 0;
random1(random1>scQtUserData.totalTrials/2) = 1;

scQtUserData.soundDeliver = random1;

%calculates delay between sound and reward
scQtUserData.soundRewDel = round(scQtUserData.soundDur-(rand(scQtUserData.totalTrials,1)*scQtUserData.waterWindow));

%calculates ITIs
k = 2.5;
p = (1-exp(-k))*rand(scQtUserData.totalTrials,1);
tau = (scQtUserData.maxITI-scQtUserData.minITI)/k;
x = round(scQtUserData.minITI + (-log(1-p))*tau)-scQtUserData.lickWindow-scQtUserData.timeDelay; 
scQtUserData.itiTime = x - scQtUserData.soundRewDel;

sendScQtControlMessage(['soundDur=',num2str(scQtUserData.soundDur)]);
sendScQtControlMessage(['rewLength',num2str(scQtUserData.rewSize)]);
sendScQtControlMessage(['lickWindow=',num2str(scQtUserData.lickWindow)]);
sendScQtControlMessage(['timeDelay=',num2str(scQtUserData.timeDelay)]);
sendScQtControlMessage(['disp(''StartSession'')']);
