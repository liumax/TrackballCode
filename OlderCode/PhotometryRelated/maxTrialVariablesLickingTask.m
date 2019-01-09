 function [trialStates, portStates, trialParams] = ...
    maxTrialVariablesLickingTask(fname)   

%These are counters
allLineNum = 0; %This is the counter for all lines. Serves as a general counter for repeats of the while loop
eventLineNum = 0;
portLineNum = 0;
trialNum = 0; %trial number counter

%These are session parameters
mouseID = [];
weight = [];
dateVal = [];
timeVal = [];
sessionID = [];
notes = [];
totalTrials = [];
toneDur = [];
toneFreq = [];
toneDB =[];
rew1Size = [];
rew2Size = [];
NASize = [];
punSize = [];
rewDelay = [];
numToneTrials = [];
numCatchTrials = [];
numFreeTrials = [];

%These are trial variables
trialType = [];
cueTimeRew1 = [];
laserTime = [];
cueTimeRew2 = [];
cueTimePun = [];
cueTimeNA = [];
playTime = [];
rewTrialTime = [];

lickTimes = [];
lickIndex = 1;

%These are variables for port-states
tStamps = [];
inStates = [];
outStates = [];


disp(fname); %This displays the fname you are using
fid = fopen(fname); %fopen opens the file
tline = fgetl(fid); %This designates tline as the next line in the logfile

while ischar(tline) %repeats loop as long as tline has characters
    findSpaces = find(tline == ' '); %looks for spaces. ~~~ lines will not have any
    if any(findSpaces) && any(str2double(tline(1:(findSpaces(1)-1)))) %activates if any spaces or any time stamp (eliminates ~~~)
        allLineNum = allLineNum + 1; %updates total line counter MAY NOT BE NECESSARY
        if isletter(tline(findSpaces(1)+1)) %searches for text after time with letters (not just portstates)
            thisTime = str2double(tline(1:(findSpaces(1)-1))); %saves the timestamp
            eventLineNum = eventLineNum + 1; %updates the eventline counter
            eventStrings{eventLineNum} = tline((findSpaces(1)+1):end); %converts the string to component of cell array
            
            %This updates trialNum counter
            if ~isempty(strfind(eventStrings{eventLineNum},'Trial = '))
                trialNum = trialNum + 1;
                
                trialType(trialNum) = 0;
                cueTimeRew1(trialNum) = 0;
                laserTime(trialNum) = 0;
                cueTimeRew2(trialNum) = 0;
                cueTimePun(trialNum) = 0;
                cueTimeNA(trialNum) = 0;
                playTime(trialNum) = 0;
                rewTrialTime(trialNum) = 0;
            end
% %             
% % %             This updates trialType so I can easily access it.
% %             if ~isempty(strfind(eventStrings{eventLineNum},'Initiating Neutral Trial'))
% %                 trialType(trialNum) = 1;
% %             elseif ~isempty(strfind(eventStrings{eventLineNum},'Initiating Reward Trial'))
% %                 trialType(trialNum) = 2;
% %             elseif ~isempty(strfind(eventStrings{eventLineNum},'Initiating Punishment Trial'))
% %                 trialType(trialNum) = 3;
% %             end

            %pull licks!
            if ~isempty(strfind(eventStrings{eventLineNum},'Lick Detected'))
                lickTimes(lickIndex,1) = str2double(tline(1:(findSpaces(1)-1)));
                lickTimes(lickIndex,2) = trialNum;
                lickIndex = lickIndex + 1;
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'Play'))
                playTime(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'RewTrial')) & isempty(strfind(eventStrings{eventLineNum},'freeRewTrials'))
                rewTrialTime(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            
            
            %The following three are for important trial landmarks,
            %including sound on, sound off, and reward delivery
            if ~isempty(strfind(eventStrings{eventLineNum},'Tone Delivered'))
                cueTimeRew1(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'LaserTrial'))
                laserTime(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            %All of the below are to populate what will basically become
            %session parameters.
            if ~isempty(strfind(eventStrings{eventLineNum},'Mouse ID'))
                mouseID = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'toneTrials'))
                numToneTrials = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'freeRewTrials'))
                numFreeTrials = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'catchTrials'))
                numCatchTrials = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'weight'))
                weight = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'date'))
                dateVal = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'time'))
                timeVal = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'sessionID'))
                sessionID = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'notes'))
                notes = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'bigReward'))
                rew1Size = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'soundDur'))
                toneDur = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'soundAmp'))
                toneDB = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'bigTone'))
                toneFreq = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'totalTrials'))
                totalTrials = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'rewDelay'))
                rewDelay = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            %to determine the type of trial. 
            if ~isempty(strfind(eventStrings{eventLineNum},'PlaySmall'))
                trialType(trialNum) = 1;
            end
            if ~isempty(strfind(eventStrings{eventLineNum},'PlayBig'))
                trialType(trialNum) = 2;
            end
            if ~isempty(strfind(eventStrings{eventLineNum},'PlayFreeRew'))
                trialType(trialNum) = 3;
            end
            if ~isempty(strfind(eventStrings{eventLineNum},'PlayBigCatch'))
                trialType(trialNum) = 4;
            end
            
            %Below code is for portstate changes
        elseif ~isnan(str2double(tline(findSpaces(1)+1))) %this picks up the portstate changes
            portLineNum = portLineNum + 1; %updates port line counter
            tStamps(portLineNum) = str2double(tline(1:(findSpaces(1)-1))); %generates tstamps record.
            inPortStateStr = ... %this line converts in states from string to double to binary
                dec2bin(str2double(tline((findSpaces(1)+1):(findSpaces(2)-1))),8);
            outPortStateStr = ...%this line converts out states from string to double to binary
                dec2bin(str2double(tline((findSpaces(2)+1):end)),8);
            inStates(portLineNum,:)=zeros(1,8); %generates space in in/out array
            outStates(portLineNum,:)=zeros(1,8);
            %This next for loop is in backwards because this is format of
            %binary output. 
            for portNum = 1:8 %fills in instates values with values from binary output. 
                inStates(portLineNum, 9-portNum) = ...
                    str2double(inPortStateStr(portNum));
                outStates(portLineNum,9-portNum) = ...
                    str2double(outPortStateStr(portNum));
            end            
        end
    end
    
    tline = fgetl(fid); %This calls the next line, in preparation for the next loop.
end

trialStates.trialType = trialType;
trialStates.cueTimeRew1 = cueTimeRew1;
trialStates.laserTime = laserTime;
trialStates.cueTimeRew2 = cueTimeRew2;
trialStates.cueTimePun = cueTimePun;
trialStates.cueTimeNA = cueTimeNA;
trialStates.playTime = playTime;
trialStates.rewTrialTime = rewTrialTime;

trialParams.mouseID = mouseID;
trialParams.weight = weight;
trialParams.dateVal = dateVal;
trialParams.timeVal = timeVal;
trialParams.sessionID = sessionID;
trialParams.notes = notes;
trialParams.rew1Size = rew1Size;
trialParams.punSize = punSize;
trialParams.toneDur  = toneDur;
trialParams.toneFreq = toneFreq;
trialParams.toneDB = toneDB;
trialParams.totalTrials = totalTrials;
trialParams.rewDelay = rewDelay;
trialParams.numToneTrials = numToneTrials;
trialParams.numCatchTrials = numCatchTrials;
trialParams.numFreeTrials = numFreeTrials;


trialParams.licking = lickTimes;

portStates.tStamps = tStamps;
portStates.inStates = inStates;
portStates.outStates = outStates;

fclose(fid);


 end
