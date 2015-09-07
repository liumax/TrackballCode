%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  
    'Session (1 = rew, 2 = rew + neut, 3 = rew + pun 4 = all)',...
    'Reward (msec):',...       
    'Punishment (msec)',...
    'Number of Trials:',...           
    'Weight:',...           
    'sessionID:',...        
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'yangDan1:';
num_lines=1;
def={'','','50','50','400','','1',''};
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
scQtUserData.session = str2num(answer{i});i=i+1;
scQtUserData.rewSize = str2num(answer{i});i=i+1;
scQtUserData.punSize = str2num(answer{i});i=i+1;
scQtUserData.totalTrials = str2num(answer{i});i=i+1;
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.taskID = 'yangDan';
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;


%hardcoded values
scQtUserData.minITI = 3000;
scQtUserData.maxITI = 6000;
scQtUserData.warning = 500;
scQtUserData.warningDelay = 1000;
scQtUserData.cueDur = 4000;
scQtUserData.graceDur = 2000;



%Placeholders
scQtUserData.lickHolder = 1;

%Date/Time
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
scQtUserData.tripSwitch = 0; %triggers end of trials at set number of trials
scQtUserData.binSize = 100; %msec per bin for rasters/histograms
scQtUserData.graphAxes = [-2:scQtUserData.binSize/1000:4.9];

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''rewSize:', num2str(scQtUserData.rewSize),''')']);
sendScQtControlMessage(['disp(''punSize:', num2str(scQtUserData.punSize),''')']);
sendScQtControlMessage(['disp(''totalTrials:', num2str(scQtUserData.totalTrials),''')']);
sendScQtControlMessage(['disp(''minITI:', num2str(scQtUserData.minITI),''')']);
sendScQtControlMessage(['disp(''maxITI:', num2str(scQtUserData.maxITI),''')']);
sendScQtControlMessage(['disp(''warning:', num2str(scQtUserData.warning),''')']);
sendScQtControlMessage(['disp(''warningDelay:', num2str(scQtUserData.cueDur),''')']);
sendScQtControlMessage(['disp(''graceDur:', num2str(scQtUserData.graceDur),''')']);
sendScQtControlMessage(['disp(''cueDur:', num2str(scQtUserData.warningDelay),''')']);
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

%calculates type of trials

scQtUserData.sessionType = zeros(scQtUserData.totalTrials,1);
if scQtUserData.session == 1
    scQtUserData.sessionType(:,1) = 2;
elseif scQtUserData.session == 2
    x = randperm(scQtUserData.totalTrials);
    y = scQtUserData.totalTrials/2;
    x(x<=y) = 1;
    x(x>y) = 2;
    scQtUserData.sessionType(:,1) = x;
elseif scQtUserData.session == 3
    x = randperm(scQtUserData.totalTrials);
    y = scQtUserData.totalTrials/2;
    x(x<=y) = 2;
    x(x>y) = 3;
    scQtUserData.sessionType(:,1) = x;
elseif scQtUserData.session == 4
    x = randperm(scQtUserData.totalTrials);
    y = scQtUserData.totalTrials/3;
    z = scQtUserData.totalTrials*2/3;
    x(x<=y) = 1;
    x(x>y & x<z)= 2;
    x(x>z) = 3;
    scQtUserData.sessionType(:,1) = x;
end

%These are temporary placeholders, which get wiped every trial
scQtUserData.LickTime = zeros (1000,1);

%These are storage for rasters
scQtUserData.neutRast = zeros (10000,2);
scQtUserData.rewRast = zeros (10000,2);
scQtUserData.punRast = zeros (10000,2);

%These are placekeeper for rasters
scQtUserData.neutHolder = 1;
scQtUserData.rewHolder = 1;
scQtUserData.punHolder = 1;

%These are storage for histograms
scQtUserData.neutHist = zeros ((scQtUserData.warningDelay + scQtUserData.cueDur + 2000)/scQtUserData.binSize,1);
scQtUserData.rewHist = zeros ((scQtUserData.warningDelay + scQtUserData.cueDur + 2000)/scQtUserData.binSize,1);
scQtUserData.punHist = zeros ((scQtUserData.warningDelay + scQtUserData.cueDur + 2000)/scQtUserData.binSize,1);

%These are holding event times for plotting
scQtUserData.cueTime = zeros(scQtUserData.totalTrials,1);
scQtUserData.rewTime = zeros(scQtUserData.totalTrials,2);
scQtUserData.punTime = zeros(scQtUserData.totalTrials,2);

sendScQtControlMessage(['rewDur=',num2str(scQtUserData.rewSize)]);
sendScQtControlMessage(['punDur=',num2str(scQtUserData.punSize)]);
sendScQtControlMessage(['graceDur=',num2str(scQtUserData.graceDur)]);
sendScQtControlMessage(['warning=',num2str(scQtUserData.warning)]);
sendScQtControlMessage(['warningDelay=',num2str(scQtUserData.warningDelay)]);
sendScQtControlMessage(['cueDur=',num2str(scQtUserData.cueDur)]);
sendScQtControlMessage(['disp(''StartSession'')']);
