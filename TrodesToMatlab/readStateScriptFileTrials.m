function trials = readStateScriptFileTrials(fileName,varargin)

%readStateScriptFileTrials(fileName,varargin)
%readStateScriptFileTrials('an123_121214_scLog.txt','inputWatch',[1:3],'trialStartTerms',{'initiate'},'trialEndTerms',{'initiate'},'checkMessageTerms',{'LEDStim'})

trialStartTerms = []; %a cell array of strings that indicate the start of a trial
trialEndTerms = []; %a cell array of strings that indicate the end of a trial
checkMessageTerms = [];
inputWatch = []; %a vector declaring which input ports to pay attention to (empty = all)
trials = {[]};

for option = 1:2:length(varargin)-1
    switch varargin{option}
       
        case 'trialStartTerms'
            trialStartTerms = varargin{option+1};
        case 'trialEndTerms'
            trialEndTerms = varargin{option+1};            
        case 'inputWatch'
            inputWatch = varargin{option+1};
        case 'checkMessageTerms'
            checkMessageTerms = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if isempty(inputWatch)
    events = parseStateScriptFile(fileName);
else 
    events = parseStateScriptFile(fileName,inputWatch);
end

inTrial = 0;
trialEventInds = [0 0];
currentTrial = 0;

%find the start and end events associated with each trial,
%based on the given strings.  Search the 'message' field for
%the strings. If the start and end terms contain the same string,
%that works too, as trials will first end, then begin on the
%same event.
for i = 1:length(events)
    if ~isempty(events(i).message)
        foundStartTerm = 0;
        foundEndTerm = 0;
        %search for end terms
        for j = 1:length(trialEndTerms)
            if ~isempty(strfind(events(i).message,trialEndTerms{j}))
                foundEndTerm = 1;               
            end
        end
        %if an end term was found, mark the event index
        if (foundEndTerm)
            if (inTrial)
                trialEventInds(currentTrial,2) = i;
            end
            inTrial = 0;           
        end
        %search for start terms
        for j = 1:length(trialStartTerms)
            if ~isempty(strfind(events(i).message,trialStartTerms{j}))
                foundStartTerm = 1;               
            end
        end
        %if a start term was found, mark the event index
        if (foundStartTerm)            
            inTrial = 1;
            currentTrial = currentTrial+1;
            trialEventInds(currentTrial,1) = i;
        end
    end
end

currentTrial = 0;
for i = 1:size(trialEventInds,1)
    currentTrial = currentTrial+1;
    currentActiveInput = 0;
    inputChangeFirstOfEach = [];
    inputChange = [];      
    outputChange = [];
    if (trialEventInds(i,2) > 0) %make sure the trial is completed
        for j = trialEventInds(i,1):trialEventInds(i,2)
            changes = find(abs(events(j).inStateDiff));
            if ((length(changes) == 1) && (changes ~= currentActiveInput))
                currentActiveInput = changes;
                inputChangeFirstOfEach = [inputChangeFirstOfEach; [events(j).time currentActiveInput]];
            elseif (length(changes) > 1)
                disp(['Warning-- more than one input changed at the same time!  Ignoring event ', num2str(i)]);
            end 
            if ~isempty(changes)
                inputChange = [inputChange;[events(j).time events(j).inStateDiff]];
            end
        end
        %messageTermsFound = zeros(1,length(checkMessageTerms));
        messageTermsFound = cell(1,length(checkMessageTerms));
        for j = trialEventInds(i,1):trialEventInds(i,2)
            changes = find(abs(events(j).outStateDiff));
            if ~isempty(changes)
                outputChange = [outputChange;[events(j).time events(j).outStateDiff]];
            end
            %check if any of the messageTerms were included in the trial
            if ~isempty(events(j).message)
                for k = 1:length(checkMessageTerms)
                    if ~isempty(strfind(events(j).message,checkMessageTerms{k}))
                        messageTermsFound{k} = [messageTermsFound{k} events(j).time];
                    end
                end
            end                 
        end
        trials{i}.timeRange = [events(trialEventInds(i,1)).time events(trialEventInds(i,2)).time]/1000;
        trials{i}.initialInputState = events(trialEventInds(i,1)).inState;
        trials{i}.initialOutputState = events(trialEventInds(i,1)).outState;
        trials{i}.inputChangeFirstOfEach = inputChangeFirstOfEach;
        trials{i}.inputChange = inputChange;
        trials{i}.outputChange = outputChange;
        %save a list of any outputs that either went high or low during the trial
        trials{i}.inputSequence = inputChangeFirstOfEach(:,2)';
        trials{i}.outputsThatChanged = find(sum(abs(outputChange(:,2:end))))';
        for j = 1:length(checkMessageTerms)
            trials{i} = setfield(trials{i},checkMessageTerms{j},messageTermsFound{j});
        end
    end    
end






