%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  
    'Reward (msec):',...       
    'Tone (3=lo,4=hi)',...
    'Number of Trials (multiple of 10):',...           
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
scQtUserData.cueDur = 500;%duration of cue
scQtUserData.rewWin = 2000; %how long the animal has to lick in response to cue
scQtUserData.histLim = [-2,5];%edges of the histogram
scQtUserData.noLick = 3000;%period where licks are not permitted
scQtUserData.percentReward = 0.45; %proportion of trials with rewarded cue delivery
scQtUserData.percentUnreward = 0.45; %proportion of trials with unrewarded cue delivery
scQtUserData.percentFree = 0.1; %proportion of trials with uncued reward

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
scQtUserData.distID = 4 - scQtUserData.cueID + 3;


pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''rewSize:', num2str(scQtUserData.rewSize),''')']);
sendScQtControlMessage(['disp(''cueID:', num2str(scQtUserData.cueID),''')']);
sendScQtControlMessage(['disp(''cueDur:', num2str(scQtUserData.cueDur),''')']);
sendScQtControlMessage(['disp(''rewWin:', num2str(scQtUserData.rewWin),''')']);
sendScQtControlMessage(['disp(''noLick:', num2str(scQtUserData.noLick),''')']);
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

%calculates which trials are which
y = randperm(scQtUserData.totalTrials);
yesRew = scQtUserData.percentReward*scQtUserData.totalTrials;
noRew = scQtUserData.percentUnreward*scQtUserData.totalTrials;
freeRew = scQtUserData.percentFree*scQtUserData.totalTrials;
y(y < yesRew) = 3;
y(y >= yesRew & y <= yesRew + noRew) = 4;
y(y > yesRew+noRew) = 5;
scQtUserData.trialDesig = y;

%These are temporary placeholders, which get wiped every trial
scQtUserData.LickTime = zeros (1000,1);

%These are storage for rasters
scQtUserData.rewRast = zeros (10000,2);
scQtUserData.distRast = zeros (10000,2);
scQtUserData.freeRast = zeros (10000,2);

%These are placekeeper for rasters
scQtUserData.rewHolder = 1;
scQtUserData.distHolder = 1;
scQtUserData.freeHolder = 1;

%These are storage for histograms
scQtUserData.rewHist = zeros(length(scQtUserData.graphAxes),1);
scQtUserData.distHist = zeros(length(scQtUserData.graphAxes),1);
scQtUserData.freeHist = zeros(length(scQtUserData.graphAxes),1);

%These are holding event times for plotting
scQtUserData.cueTime = zeros(scQtUserData.totalTrials,1); %all cue deliveries
scQtUserData.rewCueTime = zeros(scQtUserData.totalTrials,1); %reward cue deliveries
scQtUserData.distCueTime = zeros(scQtUserData.totalTrials,1); %distractor cue deliveries

%These are for first licks of specific trial types
scQtUserData.firstRew = zeros(scQtUserData.totalTrials,1);
scQtUserData.firstDist = zeros(scQtUserData.totalTrials,1);
scQtUserData.firstFree = zeros(scQtUserData.totalTrials,1);

scQtUserData.rewTime = zeros(scQtUserData.totalTrials,1); %this should plot time from cue onset to reward delivery
scQtUserData.itiRecord = zeros(scQtUserData.totalTrials,1); %this tells ITI by cue onset timing
scQtUserData.rewCounter = 0; %this counts number of rewards delivered

sendScQtControlMessage(['rewDur=',num2str(scQtUserData.rewSize)]);
sendScQtControlMessage(['cueDur=',num2str(scQtUserData.cueDur)]);
sendScQtControlMessage(['rewWin=',num2str(scQtUserData.rewWin)]);
sendScQtControlMessage(['cueID=',num2str(scQtUserData.cueID)]);
sendScQtControlMessage(['distID=',num2str(scQtUserData.distID)]);
sendScQtControlMessage(['noLick=',num2str(scQtUserData.noLick)]);
sendScQtControlMessage(['disp(''StartSession'')']);
