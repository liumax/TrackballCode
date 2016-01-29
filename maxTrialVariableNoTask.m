 function [trialStates, portStates, trialParams,trackStates] = ...
    maxTrialVariablesTuning(fname)   

%These are counters
allLineNum = 0; %This is the counter for all lines. Serves as a general counter for repeats of the while loop
eventLineNum = 0;
portLineNum = 0;
trialNum = 0; %trial number counter


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
        if ~isnan(str2double(tline(findSpaces(1)+1))) %this picks up the portstate changes
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

portStates.tStamps = tStamps;
portStates.inStates = inStates;
portStates.outStates = outStates;

trackStates.times = upAtimes;
trackStates.velocity = upA;

fclose(fid);


 end
