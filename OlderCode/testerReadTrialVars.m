 function [trialStates, portStates, trialParams] = ...
    maxTrialVariables(fname)   

%These are counters
allLineNum = 0; %This is the counter for all lines. Serves as a general counter for repeats of the while loop
eventLineNum = 0;
portLineNum = 0;
trialNum = 0; %trial number counter

%These are session parameters
fileName = [];
toneDuration = [];
toneRepetitions = [];
minITI = [];
maxITI = [];
freqStart = [];
freqEnd = [];
freqJump = [];
dbStart = [];
dbEnd = [];
dbJump = [];
totalTrials = [];

%These are trial variables
trialFreq = [];
trialDB = [];
playSound = [];
triggerMatlab = [];

%These are variables for port-states
tStamps = [];
inStates = [];
outStates = [];

%These are running speed variables

disp(fname); %This displays the fname you are using
fid = fopen(fname); %fopen opens the file
tline = fgetl(fid); %This designates tline as the next line in the logfile

while ischar(tline) %repeats loop as long as tline has characters
    findSpaces = find(tline == ' '); %looks for spaces. ~~~ lines will not have any
    findColons = find(tline == ':');
    if any(findSpaces) && any(str2double(tline(1:(findSpaces(1)-1)))) %activates if any spaces or any time stamp (eliminates ~~~)
        allLineNum = allLineNum + 1; %updates total line counter MAY NOT BE NECESSARY
        if isletter(tline(findSpaces(1)+1)) %searches for text after time with letters (not just portstates)
            thisTime = str2double(tline(1:(findSpaces(1)-1))); %saves the timestamp
            eventLineNum = eventLineNum + 1; %updates the eventline counter
            eventStrings{eventLineNum} = tline((findSpaces(1)+1):end); %converts the string to component of cell array

            %This updates trialNum counter
            if ~isempty(strfind(eventStrings{eventLineNum},'Trial = '))
                trialNum = trialNum + 1;
                trialFreq(trialNum) = 0;
                trialDB(trialNum) = 0;
                playSound(trialNum) = 0;
                triggerMatlab(trialNum) = 0;
            end

            
            %The following three are for important trial landmarks,
            %including sound on, sound off, and reward delivery
            if ~isempty(strfind(eventStrings{eventLineNum},'Play Sound'))
                playSound(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'TriggerMatlab'))
                triggerMatlab(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'Trial: '))
                trialFreq(trialNum) = str2double(tline((findColons(2)+1):findSpaces(4)-1));
                trialDB(trialNum) = str2double(tline((findColons(3)+1):end));
            end
            
            %All of the below are to populate what will basically become
            %session parameters.
            if ~isempty(strfind(eventStrings{eventLineNum},'Filename'))
                fileName = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'ToneDuration'))
                toneDuration = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'ToneRepetitions'))
                toneRepetitions = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'MinITI'))
                minITI = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'MaxITI'))
                maxITI = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'StartFrequency'))
                freqStart = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'EndFrequency'))
                freqEnd = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'OctaveJump'))
                freqJump = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'StartDB'))
                dbStart = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'EndDB'))
                dbEnd = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'dbSteps'))
                dbJump = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'TotalTrials'))
                totalTrials = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
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

trialParams.fileName = fileName;
trialParams.toneDuration = toneDuration;
trialParams.toneRepetitions = toneRepetitions;
trialParams.minITI = minITI;
trialParams.maxITI = maxITI;
trialParams.freqStart = freqStart;
trialParams.freqEnd = freqEnd;
trialParams.freqJump = freqJump;
trialParams.dbStart = dbStart;
trialParams.dbEnd = dbEnd;
trialParams.dbJump = dbJump;
trialParams.totalTrials = totalTrials;

trialStates.rewLength = trialFreq;
trialStates.preLick = trialDB;
trialStates.playSound = playSound;
trialStates.triggerMatlab = triggerMatlab;

portStates.tStamps = tStamps;
portStates.inStates = inStates;
portStates.outStates = outStates;

fclose(fid);


 end
