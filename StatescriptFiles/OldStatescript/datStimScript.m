%This is the Matlab Script

%I HAVE A HARD OFFSET FOR TIMING. This is due to a roughly 200-350ms delay
%between execution of sound(s,sf) and actual delivery of sound, as detected
%by TTL pulse. 
hardOffset = 350; %hard offset in ms.

global scQtUserData;
          
% UI prompt:
prompt = {'File Name:',...  
    'Pulse Duration (ms):',...
    'Pulse ITI (ms):',...       
    'Train ITI (min, ms):',...        
    'Train ITI (max, ms):',...
    'Pulse Number:',...
    'Train Number:',...
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'basicTuningCurve:';
num_lines=1;
def={'','10','30','8000','15000','10','300',''};
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
pulseDur = str2num(answer{i});i=i+1;
pulseITI = str2num(answer{i});i=i+1;
postPauseMin = str2num(answer{i});i=i+1;
postPauseMax = str2num(answer{i});i=i+1;
pulseNum = str2num(answer{i});i=i+1;
trainNum = str2num(answer{i});i=i+1;
scQtUserData.notes = answer{i};i=i+1;
scQtUserData.taskID = 'basicDatStim';

%Date/Time
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
scQtUserData.tripSwitch = 0; %thing to shut things down

pause(0.2);
%store total number of trials
scQtUserData.trainNum = trainNum;

%calculates ITIs
k = 2.5;
p = (1-exp(-k))*rand(trainNum,1);
tau = (postPauseMax-postPauseMin)/k;
x = round(postPauseMin + (-log(1-p))*tau); 
itiTime = x;

%save master to global
scQtUserData.ITIs = itiTime;

%% 
sendScQtControlMessage(['disp(''Filename:', scQtUserData.fileName,''')']);
sendScQtControlMessage(['disp(''PulseDuration:', num2str(pulseDur),''')']);
sendScQtControlMessage(['disp(''PulseITI:', num2str(pulseITI),''')']);
sendScQtControlMessage(['disp(''MinITI:', num2str(postPauseMin),''')']);
sendScQtControlMessage(['disp(''MaxITI:', num2str(postPauseMax),''')']);
sendScQtControlMessage(['disp(''PulseNumber:', num2str(pulseNum),''')']);
sendScQtControlMessage(['disp(''TrainNumber:', num2str(trainNum),''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

sendScQtControlMessage(['pulseNum = ',num2str(pulseNum)]);
sendScQtControlMessage(['pulseDur = ',num2str(pulseDur)]);
sendScQtControlMessage(['pulseITI = ',num2str(pulseITI)]);
sendScQtControlMessage(['pulseNum = ',num2str(pulseNum)]);

pause(1) %Need to put all my timings in before this stuff

sendScQtControlMessage(['disp(''StartSession'')']);
