%This is the Matlab Script

global scQtUserData;
          
%% UI prompt:
prompt = {'Mouse ID:',...  
    'Weight:',...  
    'Big Reward (msec):',...
    'Small Reward (msec):',...       
    'Punishment (msec):',...
    'Big Trials:',...          
    'Small Trials:',...          
    'Pun Trials:',...          
    'Free BigRew Trials:',...
    'BigRewCatch Trials:',...
    'Free Pun Trials:',...
    'PunCatch Trials:',...
    'Outcome Delay (msec):',...
    'Big Tone (port):',...       
    'Small Tone (port):',...       
    'Pun Tone (port):',...       
    'ITI (msec):',...
    'ITI (Longest) (msec):',...
    'sessionID:',...        
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'LickTask:';
num_lines=1;
def={'','','400','0','200','100','100','0','0','0','0','0','1300','','','','8000','15000','1',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

%% store all the shits
i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.weight = answer{i};i=i+1;
scQtUserData.bigRewDur = str2num(answer{i});i=i+1;
scQtUserData.smallRewDur = str2num(answer{i});i=i+1;
scQtUserData.punDur = str2num(answer{i});i=i+1;
scQtUserData.trialsBig = str2num(answer{i});i=i+1;
scQtUserData.trialsSmall = str2num(answer{i});i=i+1;
scQtUserData.trialsPun = str2num(answer{i});i=i+1;
scQtUserData.freeRew = str2num(answer{i});i=i+1;
scQtUserData.catchRewTrials = str2num(answer{i});i=i+1;
scQtUserData.freePun = str2num(answer{i});i=i+1;
scQtUserData.catchPunTrials = str2num(answer{i});i=i+1;
scQtUserData.outDelay = str2num(answer{i});i=i+1;
scQtUserData.bigTone = str2num(answer{i});i=i+1;
scQtUserData.smallTone = str2num(answer{i});i=i+1;
scQtUserData.punTone = str2num(answer{i});i=i+1;
scQtUserData.ITI = str2num(answer{i});i=i+1;
scQtUserData.ITIRange = str2num(answer{i});i=i+1;
scQtUserData.sessionID = answer{i};i=i+1;
scQtUserData.notes = answer{i};i=i+1;
scQtUserData.taskID = 'pavThreeToneWAV';

saveName = strcat('C:\Users\KreitzerLab\Desktop\',date,scQtUserData.mouseID,'RawFile');

%% now lets start calculations. 
%calculate ITIs, use exponential distribution

totalTrials = scQtUserData.trialsBig + scQtUserData.trialsSmall + scQtUserData.trialsPun + scQtUserData.freeRew + scQtUserData.catchRewTrials + scQtUserData.freePun + scQtUserData.catchPunTrials;

k = 2.5;
p = (1-exp(-k))*rand(totalTrials,1);
tau = (scQtUserData.ITIRange-scQtUserData.ITI)/k;
x = round(scQtUserData.ITI-6000 + (-log(1-p))*tau); 
scQtUserData.Master(:,1) = x;

%calculate lag for reward delivery
k = 2.5;
p = (1-exp(-k))*rand(totalTrials,1);
tau = (300)/k; %adjusted 170808 to try and improve behavior. 
x = round(scQtUserData.outDelay + (-log(1-p))*tau); 
scQtUserData.outDelayMatrix = x;

%find greatest common denominators
%to do this, we first pull all trial times, and use them to make a nx1
%vector
trialNums = [scQtUserData.trialsBig,scQtUserData.trialsSmall,scQtUserData.trialsPun,scQtUserData.freeRew,scQtUserData.catchRewTrials,scQtUserData.freePun,scQtUserData.catchPunTrials;];
%eliminate zero values, since this will break GCD code. Then recursively
%calculate GCD. Since GCD will process the GCD from the previous
%comparison, it will produce the overall GCD. 
trialNums(trialNums == 0) = [];
gcdStart = gcd(trialNums(1),trialNums(2));
for i = 2:length(trialNums)
    gcdStart = gcd(gcdStart,trialNums(i));
end

%check for errors! This would be GCD value of 1
if gcdStart == 1
    disp('Failure to Find GCD: Going Full Random Instead')
    gcdFailToggle = 1;
else
    disp(strcat('GCD Value Found',num2str(gcdStart)))
    gcdFailToggle = 0;
end

%calculate the effective sub-set size (this is the size of the sample we
%will be using the pseudorandomize within. For example, a 50 trial block
%within 200 trials total.

numIter = totalTrials/gcdStart;
divisor = gcdStart;

% trial vector will be 1 = low, 2 = hi, 3 = trials Punishment, 4= free
% reward, 5= catch rew, 5 = free punish 6 = catch punish
desigVect = zeros(numIter,1);
desigInd = 1;
desigVect(desigInd:desigInd + (scQtUserData.trialsBig/divisor)-1) = 1;desigInd = desigInd + (scQtUserData.trialsBig/divisor);
desigVect(desigInd:desigInd + (scQtUserData.trialsSmall/divisor)-1) = 2;desigInd = desigInd + (scQtUserData.trialsSmall/divisor);
desigVect(desigInd:desigInd + (scQtUserData.trialsPun/divisor)-1) = 3;desigInd = desigInd + (scQtUserData.trialsPun/divisor);
desigVect(desigInd:desigInd + (scQtUserData.freeRew/divisor)-1) = 4;desigInd = desigInd + (scQtUserData.freeRew/divisor);
desigVect(desigInd:desigInd + (scQtUserData.catchRewTrials/divisor)-1) = 5;desigInd = desigInd + (scQtUserData.catchRewTrials/divisor);
desigVect(desigInd:desigInd + (scQtUserData.freePun/divisor)-1) = 6;desigInd = desigInd + (scQtUserData.freePun/divisor);
desigVect(desigInd:desigInd + (scQtUserData.catchPunTrials/divisor)-1) = 7;desigInd = desigInd + (scQtUserData.catchPunTrials/divisor);



trialDesigs = zeros(totalTrials,1);
trialInd = 1;
for i = 1:divisor
    randSet = randperm(numIter);
    trialSet = desigVect(randSet);
    trialDesigs(trialInd:trialInd -1 + numIter) = trialSet;
    trialInd = trialInd + numIter;
end


scQtUserData.Master(:,2) = trialDesigs;


% trial vector will be 1 = low, 2 = hi, 3 = trials Punishment, 4= free
% reward, 5= catch rew, 5 = free punish 6 = catch punish

%determine rewSize order
scQtUserData.Master(:,3) = zeros(totalTrials,1);
scQtUserData.Master(trialDesigs == 1,3) = scQtUserData.smallRewDur;
scQtUserData.Master(trialDesigs == 2,3) = scQtUserData.bigRewDur;
scQtUserData.Master(trialDesigs == 3,3) = scQtUserData.punDur;
scQtUserData.Master(trialDesigs == 4,3) = 0;
scQtUserData.Master(trialDesigs == 5,3) = scQtUserData.bigRewDur;
scQtUserData.Master(trialDesigs == 6,3) = 0;
scQtUserData.Master(trialDesigs == 7,3) = scQtUserData.punDur;

%determine output port!  2 is the output port for sucrose, 3 is output for
%air puff
scQtUserData.Master(:,4) = ones(totalTrials,1);
scQtUserData.Master(trialDesigs == 1,4) = 2;
scQtUserData.Master(trialDesigs == 2,4) = 2;
scQtUserData.Master(trialDesigs == 5,4) = 2;
scQtUserData.Master(trialDesigs == 3,4) = 3;
scQtUserData.Master(trialDesigs == 7,4) = 3;
%since this uses the wavTrigger, i dont need to do any sound preparation! 

%determine sound port linkages!
scQtUserData.Master(:,5) = zeros(totalTrials,1);
scQtUserData.Master(trialDesigs == 1,5) = scQtUserData.smallTone;
scQtUserData.Master(trialDesigs == 2,5) = scQtUserData.bigTone;
scQtUserData.Master(trialDesigs == 3,5) = scQtUserData.punTone;
scQtUserData.Master(trialDesigs == 4,5) = 1;
scQtUserData.Master(trialDesigs == 5,5) = scQtUserData.bigTone;
scQtUserData.Master(trialDesigs == 6,5) = 1;
scQtUserData.Master(trialDesigs == 7,5) = scQtUserData.punDur;

%store information about time/date
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));
scQtUserData.tripSwitch = 0;

% my additional fields:
scQtUserData.trial = 0; % keep track of trial number
scQtUserData.failTrig = 0; %trigger for failure of the graphing code. 

pause(0.5);

%% Display these to statescript log file for storage. 
sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''weight:', scQtUserData.weight,''')']);
sendScQtControlMessage(['disp(''bigReward:', num2str(scQtUserData.bigRewDur),''')']);
sendScQtControlMessage(['disp(''smallReward:', num2str(scQtUserData.smallRewDur),''')']);
sendScQtControlMessage(['disp(''punDur:', num2str(scQtUserData.punDur),''')']);
sendScQtControlMessage(['disp(''trialsBig:', num2str(scQtUserData.trialsBig),''')']);
sendScQtControlMessage(['disp(''trialsSmall:', num2str(scQtUserData.trialsSmall),''')']);
sendScQtControlMessage(['disp(''trialsPun:', num2str(scQtUserData.trialsPun),''')']);
sendScQtControlMessage(['disp(''freeRewTrials:', num2str(scQtUserData.freeRew),''')']);
sendScQtControlMessage(['disp(''catchRewTrials:', num2str(scQtUserData.catchRewTrials),''')']);
sendScQtControlMessage(['disp(''freePunTrials:', num2str(scQtUserData.freePun),''')']);
sendScQtControlMessage(['disp(''catchPunTrials:', num2str(scQtUserData.catchPunTrials),''')']);
sendScQtControlMessage(['disp(''outDelay:', num2str(scQtUserData.outDelay),''')']);
sendScQtControlMessage(['disp(''bigTone:', num2str(scQtUserData.bigTone),''')']);
sendScQtControlMessage(['disp(''smallTone:', num2str(scQtUserData.smallTone),''')']);
sendScQtControlMessage(['disp(''punTone:', num2str(scQtUserData.punTone),''')']);
sendScQtControlMessage(['disp(''ITIShort:', num2str(scQtUserData.ITI),''')']);
sendScQtControlMessage(['disp(''ITILong:', num2str(scQtUserData.ITIRange),''')']);
sendScQtControlMessage(['disp(''taskID:', scQtUserData.taskID,''')']);
sendScQtControlMessage(['disp(''date:', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''time:', scQtUserData.time,''')']);
sendScQtControlMessage(['disp(''sessionID:', scQtUserData.sessionID,''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

pause(1) %Need to put all my timings in before this stuff

%% generate space in structure for storage of information that I care about!
%This is for soundOn times. This allows for calculations of all the other
%things!
scQtUserData.soundOn = zeros(scQtUserData.toneTrials,1);
%this is for licking latency
scQtUserData.lickLat = zeros(scQtUserData.toneTrials,1);


%variables for tracking licking. Each entry here will be a combination of
%the lick time relative to the sound, the trial number, and the type of
%trial.
scQtUserData.licks = zeros(1000,4);
scQtUserData.lickCounter = 1;

scQtUserData.lickHist = zeros(80,2); %This is optimized for looking at an 8 second window
%with 2 sec before sound onset, sound, and 3 seconds after. Set for 100 ms
%bins. First column is for small reward, second column for big rewards.

%toggle for triggering storage of information when plotting fails.
scQtUserData.PlotToggle = 0;

%save all the input data
save(saveName,'scQtUserData')


scQtUserData.lickAxes = [-2:0.1:5.9]; %axis for histogram
%send initial information to the mbed
sendScQtControlMessage(['toneOutDel =',num2str(scQtUserData.RewDelayMatrix(1))]);
% sendScQtControlMessage(['signalDel =3000']); %this is the delay after reward delivery before triggering next thing. 
sendScQtControlMessage(['disp(''StartSession'')']);
