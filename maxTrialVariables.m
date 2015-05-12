 function [trialStates, portStates, trialParams] = ...
    maxTrialVariables(fname)   

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
soundDur = [];
minRew = [];
maxRew = [];
rewProb = [];
laserProb = [];
blocks = [];
blockSize = [];
minITI = [];
maxITI = [];
waterWindow = [];
baitDur = [];
training = [];

%These are trial variables
rewLength = [];
preLick = [];
antiLick = [];
postLick = [];
otherLick = [];
soundOn = [];
soundOff = [];
rewDelivery = [];

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
                preLick(trialNum) = 0;
                antiLick(trialNum) = 0;
                postLick(trialNum) = 0;
                otherLick(trialNum) = 0;
                soundOn(trialNum) = 0;
                soundOff(trialNum) = 0;
                rewDelivery(trialNum) = 0;
            end
            
            %This updates rewlength so I can easily access it.
            if ~isempty(strfind(eventStrings{eventLineNum},'rewLength'))
                rewLength(trialNum) = str2double(eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}=='=')+1):end));
            end
            
            %The following three are for important trial landmarks,
            %including sound on, sound off, and reward delivery
            if ~isempty(strfind(eventStrings{eventLineNum},'SoundOn'))
                soundOn(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'SoundOff'))
                soundOff(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'Reward Delivered'))
                rewDelivery(trialNum) = str2double(tline(1:(findSpaces(1)-1)));
            end
            
            %This updates lick counters for all three classes of licks. 
            if ~isempty(strfind(eventStrings{eventLineNum},'trialState = 1'))
                preLick(trialNum) = preLick(trialNum) + 1;
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'trialState = 2'))
                antiLick(trialNum) = antiLick(trialNum) + 1;
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'trialState = 3'))
                postLick(trialNum) = postLick(trialNum) + 1;
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'trialState = 4'))
                otherLick(trialNum) = otherLick(trialNum) + 1;
            end
            
            %All of the below are to populate what will basically become
            %session parameters.
            if ~isempty(strfind(eventStrings{eventLineNum},'Mouse ID'))
                mouseID = eventStrings{eventLineNum}(...
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
            
            if ~isempty(strfind(eventStrings{eventLineNum},'soundDur'))
                soundDur = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'minRew'))
                minRew = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'maxRew'))
                maxRew = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'rewProb'))
                rewProb = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'laserProb'))
                laserProb = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'blocks'))
                blocks = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'blockSize'))
                blockSize = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'minITI'))
                minITI = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'maxITI'))
                maxITI = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'waterWindow'))
                waterWindow = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'baitDur'))
                baitDur = eventStrings{eventLineNum}(...
                    (find(eventStrings{eventLineNum}==':')+1):end);
            end
            
            if ~isempty(strfind(eventStrings{eventLineNum},'Training'))
                training = eventStrings{eventLineNum}(...
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

trialStates.rewLength = rewLength;
trialStates.preLick = preLick;
trialStates.antiLick = antiLick;
trialStates.postLick = postLick;
trialStates.otherLick = otherLick;
trialStates.soundOn = soundOn;
trialStates.soundOff = soundOff;
trialStates.rewDelivery = rewDelivery;

trialParams.mouseID = mouseID;
trialParams.weight = weight;
trialParams.dateVal = dateVal;
trialParams.timeVal = timeVal;
trialParams.sessionID = sessionID;
trialParams.notes = notes;
trialParams.soundDur = soundDur;
trialParams.minRew = minRew;
trialParams.maxRew = maxRew;
trialParams.rewProb = rewProb;
trialParams.laserProb = laserProb;
trialParams.blocks = blocks;
trialParams.blockSize = blockSize;
trialParams.minITI = minITI;
trialParams.maxITI = maxITI;
trialParams.waterWindow = waterWindow;
trialParams.baitDur = baitDur;
trialParams.training = training;

portStates.tStamps = tStamps;
portStates.inStates = inStates;
portStates.outStates = outStates;

fclose(fid);


 end
