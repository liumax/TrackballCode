%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  1
    'Reward (msec):',...       2
    'Trials (even #):',...           5
    'Weight:',...           8
    'sessionID:',...        9
    'Notes:'}; %the bracket is to end the prompt     10
dlg_title = 'StopTask:';
num_lines=1;
def={'','50','200','','1',''};
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
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.taskID = 'LickTask';
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;

scQtUserData.minITI=3000;
scQtUserData.maxITI=scQtUserData.minITI+6000;


%look at these later
scQtUserData.waitPeriod = 1000; 

scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.tripSwitch = 0;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
% scQtUserData.trInit = false;
% scQtUserData.trBait = false;

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''rewSize:', num2str(scQtUserData.rewSize),''')']);
sendScQtControlMessage(['disp(''totalTrials:', num2str(scQtUserData.totalTrials),''')']);
sendScQtControlMessage(['disp(''minITI:', num2str(scQtUserData.minITI),''')']);
sendScQtControlMessage(['disp(''maxITI:', num2str(scQtUserData.maxITI),''')']);
sendScQtControlMessage(['disp(''weight:', scQtUserData.weight,''')']);
sendScQtControlMessage(['disp(''taskID:', scQtUserData.taskID,''')']);
sendScQtControlMessage(['disp(''date:', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''time:', scQtUserData.time,''')']);
sendScQtControlMessage(['disp(''sessionID:', scQtUserData.sessionID,''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

pause(1) %Need to put all my timings in before this stuff

%calculates ITIs
k = 2.5;
p = (1-exp(-k))*rand(scQtUserData.totalTrials,1);
tau = (scQtUserData.maxITI-scQtUserData.minITI)/k;
x = round(scQtUserData.minITI + (-log(1-p))*tau); 
scQtUserData.itiTime = x;

%This is for lightOn times
scQtUserData.lightOn = zeros(scQtUserData.totalTrials,1);

scQtUserData.stopTime = zeros(scQtUserData.totalTrials,1);

scQtUserData.rewardTime = zeros(scQtUserData.totalTrials,1);

scQtUserData.latTime = zeros(scQtUserData.totalTrials,1);

%this is for ITIs
scQtUserData.ITI = zeros(scQtUserData.totalTrials,1);


sendScQtControlMessage(['rewLength=',num2str(scQtUserData.rewSize)]);
sendScQtControlMessage(['waitPeriod=',num2str(scQtUserData.waitPeriod)]);
sendScQtControlMessage(['disp(''StartSession'')']);
