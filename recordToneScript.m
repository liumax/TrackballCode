%This is the Matlab Script

global scQtUserData;


% UI prompt:
prompt = {'Mouse ID:',...  
    'Drive ID:',...       
    'Tone Duration (ms):',...           
    'Min ITI (s):',...
    'Max ITI (s):',...
    'Number of Recording Sessions:',...
    'Stim Pulse Duration (ms):',...
    'Stim Frequency (Hz):',...
    'Stim Number:'...
    'Stim Timing (before to cue, ms):',...
    'Frequency of Stimulation (%):'...
    'Notes:'}; %the bracket is to end the prompt     1
dlg_title = 'StopTask:';
num_lines=1;
def={'','','500','55','90','100','5','25','20','100','20',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.driveID = answer{i};i=i+1;
scQtUserData.toneDur = str2num(answer{i});i=i+1;
scQtUserData.minITI = str2num(answer{i})*1000;i=i+1;
scQtUserData.maxITI = str2num(answer{i})*1000;i=i+1;
scQtUserData.recordNum = str2num(answer{i})*1000;i=i+1;
scQtUserData.stimDur = str2num(answer{i});i=i+1;
scQtUserData.stimFreq = 1/str2num(answer{i})*1000;i=i+1;
scQtUserData.stimNum = str2num(answer{i});i=i+1;
scQtUserData.stimTime = str2num(answer{i});i=i+1;
scQtUserData.freqStim = round(1/(str2num(answer{i})/100));i=i+1;
scQtUserData.notes = answer{i};i=i+1;

scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.triggerSwitch = 0;
scQtUserData.triggerCounter = 1;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
scQtUserData.maxTrials = 10000; %this is for calculating ITIs. Placeholder of obscenely large size to avoid issues of running out of times.

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''Drive ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''Tone Duration:', num2str(scQtUserData.toneDur),''')']);
sendScQtControlMessage(['disp(''minITI:', num2str(scQtUserData.minITI),''')']);
sendScQtControlMessage(['disp(''maxITI:', num2str(scQtUserData.maxITI),''')']);
sendScQtControlMessage(['disp(''Number of Recording Sessions:', num2str(scQtUserData.recordNum),''')']);
sendScQtControlMessage(['disp(''Stim Pulse Duration:', num2str(scQtUserData.stimDur),''')']);
sendScQtControlMessage(['disp(''Stim Frequency:', num2str(scQtUserData.stimFreq),''')']);
sendScQtControlMessage(['disp(''Stim Number:', num2str(scQtUserData.stimNum),''')']);
sendScQtControlMessage(['disp(''Stim Timing:', num2str(scQtUserData.stimTime),''')']);
sendScQtControlMessage(['disp(''Percent of Trials for Stim:', num2str(scQtUserData.freqStim),''')']);
sendScQtControlMessage(['disp(''date:', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

pause(1) %Need to put all my timings in before this stuff

%calculates ITIs
k = 2.5;
p = (1-exp(-k))*rand(scQtUserData.maxTrials,1);
tau = (scQtUserData.maxITI-scQtUserData.minITI)/k;
x = round(scQtUserData.minITI + (-log(1-p))*tau); 
scQtUserData.itiTime = x;

sendScQtControlMessage(['toneDur=',num2str(scQtUserData.toneDur)]);
sendScQtControlMessage(['stimDur=',num2str(scQtUserData.stimDur)]);
sendScQtControlMessage(['stimFreq=',num2str(scQtUserData.stimFreq)]);
sendScQtControlMessage(['stimNum=',num2str(scQtUserData.stimNum)]);
sendScQtControlMessage(['stimTime=',num2str(scQtUserData.stimTime)]);
sendScQtControlMessage(['disp(''StartSession'')']);
