%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  
    'Reward (msec):',...       
    'Tone (0=lo,1=hi)',...
    'Number of Trials:',...           
    'Weight:',...           
    'sessionID:',...        
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'auditoryLickOperantTraining:';
num_lines=1;
def={'','50','','400','','1',''};
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
scQtUserData.cueID = str2num(answer{i});i=i+1;
scQtUserData.totalTrials = str2num(answer{i});i=i+1;
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.taskID = 'auditoryLickOperantTraining';
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;


%hardcoded values
scQtUserData.minITI = 7000;
scQtUserData.maxITI = 12000;
scQtUserData.cueDelay = 500;
scQtUserData.histLim = [-2,4.9];%edges of the histogram

%Placeholders
scQtUserData.lickHolder = 1;

%Date/Time
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
scQtUserData.tripSwitch = 0; %triggers end of trials at set number of trials
scQtUserData.binSize = 100; %msec per bin for rasters/histograms
scQtUserData.graphAxes = [scQtUserData.histLim(1):scQtUserData.binSize/1000:scQtUserData.histLim(2)];


pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''rewSize:', num2str(scQtUserData.rewSize),''')']);
sendScQtControlMessage(['disp(''cueID:', num2str(scQtUserData.cueID),''')']);
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


%These are temporary placeholders, which get wiped every trial
scQtUserData.LickTime = zeros (1000,1);

%These are storage for rasters

scQtUserData.rewRast = zeros (10000,2);

%These are placekeeper for rasters

scQtUserData.rewHolder = 1;

%These are storage for histograms
scQtUserData.rewHist = zeros ((7000)/scQtUserData.binSize,1);%%% PROBLEM HERE???

%These are holding event times for plotting
scQtUserData.cueTime = zeros(scQtUserData.totalTrials,1);
scQtUserData.rewTime = zeros(scQtUserData.totalTrials,1);

sendScQtControlMessage(['rewDur=',num2str(scQtUserData.rewSize)]);
sendScQtControlMessage(['cueID=',num2str(scQtUserData.cueID)]);
sendScQtControlMessage(['disp(''StartSession'')']);
