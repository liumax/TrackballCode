
%inputs: rotary1 and rotary2 are the single column vectors of the MBED
%portstates for inputs. timeStamps is the timestamps from the MBED
%portstates, converted TO SECONDS. interpStep is the interpolation step,
%which should be set to 0.1 seconds. 

function [x] = functionMBEDrotary(rotary1,rotary2,timeStamps,interpStep);
% sampRate = 30000; %trodes sampling rate.
% interpStep = 0.1; %step size in seconds, for the interpolation. 
%for current MBED (170607), 4 is rotary1, 5 is rotary2

%steps are hard coded values based on sample data that are meant to
%compensate for differences in how the rotary encoder can compute the
%distance between up to up, up to down, down to down, and down to up states
%between the two photodiodes.
stepNorm = 0.5;
stepBig = 0.572;
stepSmall = 0.316;

%% Data Extraction and Preprocessing: 
%The columns are times, the two inputs, and the state, which is a
%combination of the different inputs. 
timeStateArray = zeros(length(timeStamps),4);
timeStateArray(:,1) = timeStamps;
timeStateArray(:,2) = rotary1;
timeStateArray(:,3) = rotary2;
%the (2:end) fudge factor for dio4States is there to compensate for the
%fact taht the first datapoint is at time zero, which is lost in the unique
%function call. 

%now go through this array, and systematically fill placeholders (5s) with
%previous state. Also fill out forth column, which represents distinct
%states.
for dioInd = 1:length(timeStamps)
    
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
    timeStateArray(repFinder,:) = [];
end

catTimes = timeStateArray(:,1);
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
revRev2For = unique([strfind(stateString,'818'),strfind(stateString,'181')]);

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
newTimes = [(round(timeStateArray(1,1)*(1/interpStep)))*interpStep:interpStep:(round(timeStateArray(end,1)*(1/interpStep)))*interpStep];

%first, basically map my current information onto newTimes without
%interpolation
fillDist(1) = 0;
for intInd = 2:length(newTimes)
    %check to see if distance has changed at all
    changeFind = find(timeStateArray(:,1) >= newTimes(intInd - 1) & timeStateArray(:,1) <= newTimes(intInd));
    if isempty(changeFind)
        %if the time hasnt passed some actualy change, do nothing
        fillDist(intInd) = fillDist(intInd-1);
    else
        %if there are things that have happened, I need to fill in these
        %data points. I assume with a 10ms interpolation, any events of
        %siginificance are going to be things that can be glommed together
        %into a single event. 
        totalChange = sum(distInc(changeFind)); 
        fillDist(intInd) = fillDist(intInd-1) + totalChange;
    end
    
end

% %this is tester code. Based on these tests, I think I'm going to stick with
% %the 31 smooth window and the lowess method. 
% smoothDist(:,1) = smooth(fillDist,11);
% smoothDist(:,2) = smooth(fillDist,11,'lowess');
% smoothDist(:,3) = smooth(fillDist,31);
% smoothDist(:,4) = smooth(fillDist,31,'lowess');
% 
% testVel(:,1) = diff(smoothDist(:,1));
% testVel(:,2) = diff(smoothDist(:,2));
% testVel(:,3) = diff(smoothDist(:,3));
% testVel(:,4) = diff(smoothDist(:,4));

% Tested various interpolation methods, find that PCHIP seems to produce
% most believable result out of all possibilities. 

try
    newDist = smooth(fillDist,3,'lowess');
%     newDist = smooth(fillDist,31,'lowess'); %170707 Adjusting to 7. 31 is
%     unnecessary levels of smoothing... 170811 adjusting to 3. 
    failTrigger = 0;
catch
    disp('Distance Smoothing Failed. Suggests there is no movement. Replacing with zeros')
    cumDist = zeros(length(catTimes),1);
    newDist = zeros(length(newTimes),1);
    %also trigger failure notification for downstream code
    failTrigger = 1;
end

%insert fail trigger if overall cumulative distance is zero
if cumDist(end) == 0
    failTrigger = 1;
    disp('No Coherent Movement Detected, Do Have Wobbles')
end

mouseVel = diff(newDist)/interpStep;
velTimes = newTimes(1:end-1); %this adds a fudge factor because of the change in the length of the array.

%reshape and make distance and time vectors into two column arrays
catTimes = reshape(catTimes,[],1);
cumDist = reshape(cumDist,[],1);

newTimes = reshape(newTimes,[],1);
newDist = reshape(newDist,[],1);

velTimes = reshape(velTimes,[],1);
mouseVel = reshape(mouseVel,[],1);

%find start times of velocity

locoBinary = zeros(length(mouseVel),1);
locoBinary(abs(mouseVel)>1) = 1; 

locoStarts = find(diff(locoBinary) == 1)+1;
locoEnds = find(diff(locoBinary) == -1);

loopTrig = 0;
if failTrigger == 0
    while loopTrig == 0;
        %combine these two to generate a single timeline. 
        [fullLocoTrace,iLoco,iC] = unique([locoStarts;locoEnds],'rows');
        %since locoStarts come first, I can separate iC out into starts and
        %ends. 
        iStarts = iC(1:length(locoStarts));
        iEnds = iC(length(locoStarts)+1:end);
        %now i need to see if the animal was locomoting to being with. If so, then
        %the first even will be a locomotion end without a preceding locomotion
        %start.
        if iEnds(1) == 1
            locoStarts(2:end+1) = locoStarts(1:end);
            locoStarts(1) = 1;
        end

        %now i need to check and see if the last value is an end or a start. If a
        %start, then I need to insert an end at the final value for locomotion. 
        if iStarts(end) == length(fullLocoTrace);
            %add additional time point to the end of locoEnds.
            locoEnds(end+1) = length(locoBinary);
        end

        if iStarts(1) == 1 & iEnds(end) == length(fullLocoTrace);
            loopTrig = 1;
        end

    end
elseif failTrigger == 1
    iStarts = 0;
    iEnds = 0;
end

x=struct;
x.Distance = [newTimes,newDist];
x.Velocity = [velTimes,mouseVel];
x.RawDistance = [catTimes,cumDist];
x.RawData = timeStateArray;
x.BinaryLocomotion = locoBinary;
% x.BinaryLocoRev = revBinary;
x.LocoStarts = velTimes(locoStarts);
x.LocoEnds = velTimes(locoEnds);


end


