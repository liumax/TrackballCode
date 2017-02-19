

function [s] = functionRotaryExtraction(s,sampRate,interpStep);
sampRate = 30000; %trodes sampling rate.
interpStep = 0.01; %step size in seconds, for the interpolation. 

%steps are hard coded values based on sample data that are meant to
%compensate for differences in how the rotary encoder can compute the
%distance between up to up, up to down, down to down, and down to up states
%between the two photodiodes.
stepNorm = 0.5;
stepBig = 0.572;
stepSmall = 0.316;

%% Data Extraction and Preprocessing: 
%find the targeted files, extract file names
[D3FileName] = functionFileFinder(subFoldersCell,'DIO','D3');
D3FileName = D3FileName{1};

[D4FileName] = functionFileFinder(subFoldersCell,'DIO','D4');
D4FileName = D4FileName{1};

%pull data from targeted files.
[DIO3Data] = readTrodesExtractedDataFile(D3FileName);
[DIO4Data] = readTrodesExtractedDataFile(D4FileName);

%turns out i wired things up funny, DIO4 is actually the first one to cross
%the plastic in forward motion. XD

%pull times and states. 
dio3Times = double(DIO3Data.fields(1).data)/sampRate;
dio3States = double(DIO3Data.fields(2).data);

dio4Times = double(DIO4Data.fields(1).data)/sampRate;
dio4States = double(DIO4Data.fields(2).data);

%the first values are initialized values when recording starts. These will
%help track the direction of the wheel.
dio3Init = dio3States(1);
dio4Init = dio4States(1); %NOTE ADJUSTED FROM 2

%now, generate a concatenated form of all time points. ia and ic give
%indices, where catTimes = dioTimes(ia), and dioTimes = catTimes(ic)
[catTimes ia ic] = unique([dio3Times;dio4Times]);

timeLim = length(dio3Times); %to differentiate between dio3 and dio4 times

%here, ia has the indexing to determine how things should be
%separated/indexed. Here, lets figure out which of the things belong to
%which DIO channel
indexDio1 = find(ia<=timeLim);
indexDio2 = find(ia>timeLim);

%dio2 indices are screwed up, since they are offset by the length of
%dio1times. fix this.
trueIndDio2 = ia(indexDio2)-timeLim;

%now put everything into a big array of all time points. This should
%generate an array of length equal to all time stamps. This fills in the
%array with data from the dio states, leaving excess spaces with a
%placeholder (5). 

%The columns are times, dio3, dio4, and state, which is defined as one of
%four combinations of dio3 and dio4.
timeStateArray = zeros(length(catTimes),4);
timeStateArray(:,1) = catTimes;
timeStateArray(:,2:3) = 5; %5 serves as a placeholder to show sites that have no change
timeStateArray(1,2) = dio3Init;
timeStateArray(1,3) = dio4Init;
timeStateArray(indexDio1,2) = dio3States;
timeStateArray(indexDio2,3) = dio4States(2:end); 
%the (2:end) fudge factor for dio4States is there to compensate for the
%fact taht the first datapoint is at time zero, which is lost in the unique
%function call. 

%now go through this array, and systematically fill placeholders (5s) with
%previous state. Also fill out forth column, which represents distinct
%states.
for dioInd = 1:length(catTimes)
    %check if have a 5 value for either. If so, this must be erased and replaced with
    %last value.
    if dioInd > 1
        if timeStateArray(dioInd,2) == 5
            timeStateArray(dioInd,2) = timeStateArray(dioInd-1,2);
        else
            timeStateArray(dioInd,3) = timeStateArray(dioInd-1,3);
        end
    end
    
    %now need logic to fill out states!
    if timeStateArray(dioInd,2) == 1 & timeStateArray(dioInd,3) == 1
        timeStateArray(dioInd,4) = 1;
    elseif timeStateArray(dioInd,2) == 0 & timeStateArray(dioInd,3) == 0
        timeStateArray(dioInd,4) = 2;
    elseif timeStateArray(dioInd,2) == 0 & timeStateArray(dioInd,3) == 1
        timeStateArray(dioInd,4) = 4;
    elseif timeStateArray(dioInd,2) == 1 & timeStateArray(dioInd,3) == 0
        timeStateArray(dioInd,4) = 8;
    end
end

%eliminate entire row if there are duplicate values of the same state
repFinder = find(diff(timeStateArray(:,4))==0);
if ~isempty(repFinder)
    timeStateArray(4,:) = [];
end
%% Process State Data into Distance
%process state data by converting states into a string, which allows for
%easy search via strfind. 

%convert states into string
stateString = num2str(timeStateArray(:,4)');
%delete spaces in string
stateString(stateString == ' ') = [];

%find all cases where running forward
runFor = unique([strfind(stateString,'4182'),...
    strfind(stateString,'1824'),...
    strfind(stateString,'8241'),...
    strfind(stateString,'2418')]);
%find all cases of clear running backwards
runRev = unique([strfind(stateString,'8142'),...
    strfind(stateString,'1428'),...
    strfind(stateString,'4281'),...
    strfind(stateString,'2814')]);
%find all cases of wobble
wobble1 = unique([strfind(stateString,'828'),...
    strfind(stateString,'282')]);
wobble2 = unique([strfind(stateString,'424'),...
    strfind(stateString,'242')]);
%find all cases of reversal
revFor2Rev = unique([strfind(stateString,'414'),strfind(stateString,'141')]);
revRev2For = unique([strfind(stateString,'818'),strfind(stateString,'818')]);

%store these data in an array of size equal to the number of state changes
%by the types of states available
bigArray = zeros(length(stateString),6);

for dioInd = 1:length(runFor)
    bigArray(runFor(dioInd):runFor(dioInd)+3,1) = 1;
end

for dioInd = 1:length(runRev)
    bigArray(runRev(dioInd):runRev(dioInd)+3,2) = 1;
end

for dioInd = 1:length(wobble1)
    bigArray(wobble1(dioInd):wobble1(dioInd)+2,3) = 1;
end

for dioInd = 1:length(wobble2)
    bigArray(wobble2(dioInd):wobble2(dioInd)+2,4) = 1;
end

for dioInd = 1:length(revFor2Rev)
    bigArray(revFor2Rev(dioInd):revFor2Rev(dioInd)+2,5) = 1;
end

for dioInd = 1:length(revRev2For)
    bigArray(revRev2For(dioInd):revRev2For(dioInd)+2,6) = 1;
end

%confirm that the big array actually accounts for all events. 
tester = sum(bigArray,2);
tester = find(tester == 0);
if length(tester) == 1
    disp('One Event Is Unaccounted for, Proceeding')
elseif isempty(tester)
    disp('All Points Accounted For, Proceeding')
else
    disp('High Number of Time Points Unaccounted For')
end

%delete stored variables.
runFor = [];
runRev = [];
wobble1 = [];
wobble2 = [];
revFor2Rev = [];
revRev2For = [];
tester = [];

%insert distances. Only count coherent movement in either direction
distInc = zeros(length(bigArray),1);
distInc(bigArray(:,1) == 1) = stepNorm;
distInc(bigArray(:,2) == 1) = -stepNorm;
%replace values for non-0.5 increments
distInc(bigArray(:,1) == 1 & timeStateArray(:,4) == 2) = stepSmall;
distInc(bigArray(:,1) == 1 & timeStateArray(:,4) == 1) = stepBig;
distInc(bigArray(:,2) == 1 & timeStateArray(:,4) == 2) = -stepSmall;
distInc(bigArray(:,2) == 1 & timeStateArray(:,4) == 1) = -stepBig;

%calculate cumulative sum for distance traveled
cumDist = cumsum(distInc);

%What we want to do now is interpolate onto a new time frame, so that we
%can fill in the spaces between actual data entries. 

%generate new time frame
newTimes = [(round(catTimes(1)*(1/interpStep)))*interpStep:interpStep:(round(catTimes(end)*(1/interpStep)))*interpStep];

% Tested various interpolation methods, find that PCHIP seems to produce
% most believable result out of all possibilities. 

newDist = interp1(catTimes,cumDist,newTimes,'pchip');

mouseVel = diff(newDist);
velTimes = newTimes(1:end-1); %this adds a fudge factor because of the change in the length of the array.

%reshape and make distance and time vectors into two column arrays
catTimes = reshape(catTimes,[],1);
cumDist = reshape(cumDist,[],1);

newTimes = reshape(newTimes,[],1);
newDist = reshape(newDist,[],1);

velTimes = reshape(velTimes,[],1);
mouseVel = reshape(mouseVel,[],1);


x=struct;
x.Distance = [newTimes,newDist];
x.Velocity = [velTimes,mouseVel];
x.RawDistance = [catTimes,cumDist];
x.RawData = timeStateArray;

s.RotaryData = x;

end


