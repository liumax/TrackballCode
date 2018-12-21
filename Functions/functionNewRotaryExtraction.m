%processes rotary encoder code. 
%Inputs: name1 (name for Din7 for Tony, Din3 for Max)
% name2 (name for Din8 for Tony, Din4 for Max);

%Outputs: Outputs a structured array with the following: 
    %Distance: two column vector with time points in trodes samples and
    %distance in cm
    %Velocity: velocity in cm/sec. 

function [funcOut] = functionNewRotaryExtraction(name1,name2);
% sampRate = 30000; %trodes sampling rate.
% interpStep = 0.001; %step size in seconds, for the interpolation. 

%Since things are encoded in quadrature, the state transition represents
%the movement of 1/800 of the circumference. Radius for disk is
%approximately 7cm, so diameter is 14. This means radius is 44 cm. 44 cm
%divided into 800ths is 0.055 cm.
stepSize = 0.055;
% stepBig = 0.3;
% stepSmall = 0.3;

%% Data Extraction and Preprocessing: 
disp('Extracting Data from DIO Files')
%pull data from targeted files.
[DIO3Data] = readTrodesExtractedDataFile(name1);
[DIO4Data] = readTrodesExtractedDataFile(name2);
disp('DIO Data Extracted')
%pull times and states. 
dio3Times = double(DIO3Data.fields(1).data);
dio3States = double(DIO3Data.fields(2).data);

dio4Times = double(DIO4Data.fields(1).data);
dio4States = double(DIO4Data.fields(2).data);

timeMin = min([dio3Times;dio4Times])-10;
timeMax = max([dio3Times;dio4Times])+ 10;

%the first values are initialized values when recording starts. These will
%help track the direction of the wheel.
dio3Init = dio3States(1);
dio4Init = dio4States(1); %NOTE ADJUSTED FROM 2

%170726 Having problem where there is probably overlap of times between
%dio3 and dio4. Trying to resolve with code that detects the same time
%across both and deletes extraneous pulses.
% timeMin = dio3Times(1);
% timeMax = dio3Times(end);
[C,ia,ib] = intersect(dio3Times,dio4Times);
if length(ia) > 1 | length(ib) >1
    %delete first value, because this is always the first time point
    ia(1) = [];
    ib(1) = [];
    dio3Times(ia) = [];
    dio3States(ia) = [];
    dio4Times(ia) = [];
    dio4States(ia) = [];
end

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
disp('State Array Generated')

%eliminate entire row if there are duplicate values of the same state
repFinder = find(diff(timeStateArray(:,4))==0);
if ~isempty(repFinder)
    timeStateArray(repFinder,:) = [];
end

%update catTimes
catTimes = timeStateArray(:,1);
%% Process State Data into Distance

% stateArray = [0,0,-1,1;
%     0,0,1,-1;
%     1,-1,0,0;
%     -1,1,0,0];
%because the difference of each state produces a unique value, we can
%assign different changes to those values. 

%generate changes in state
timeStateDiff = diff(timeStateArray(:,4));

distChange = zeros(length(timeStateDiff),1);

%find all positive steps forward
negSteps= find(timeStateDiff == -7 | timeStateDiff == -2 | timeStateDiff == 3 | timeStateDiff == 6);
posSteps = find(timeStateDiff == -3 | timeStateDiff == -6 | timeStateDiff == 2 | timeStateDiff == 7);
distChange(posSteps)  =  stepSize;
distChange(negSteps)  =  -stepSize;

disp('Replaced States with Distance Changes')
%fill in time points of there is a mismatch of DIO times. Insert extra data
%points!

if timeMin < catTimes(1)
    %shift catTimes by 1. 
    catTimes(2:end+1) = catTimes(1:end);
    %fill in first point with minimum time
    catTimes(1) = timeMin;
    distChange(2:end+1) = distChange(1:end);
    distChange(1) = 0;

end

if timeMax > catTimes(end)
    %insert extra value for catTimes and cumDist
    catTimes(end+1) = timeMax;
    distChange(end+1) = 0;
end

% Fill in extra time points. We want this to basically fill in time points 
%so that any points after a registered state change reflect the previous 
%state. This way, we can fill in spaces between long gaps with no movement.
newTimes = catTimes - catTimes(1) + 1;
newVel = zeros(timeMax-timeMin+1,1);
for i = 1:length(distChange)
    newVel(newTimes(i+1)) = distChange(i); %note that this uses the next time point. This may not be the best method. Could use preceding time point.
end

%now convert to distance
newDist = cumsum(newVel);
%smooth over 50 ms window
smoothDist = smooth(newDist,1501);
%now downsample to single ms resolution.
smoothDist = downsample(smoothDist,30);
if smoothDist(end) < 100
    failTrigger = 1;
else
    failTrigger = 0;
end
%generate appropriate time vector. 
smoothDistTime = [timeMin:30:timeMin+30*(length(smoothDist)-1)];

%we now have a distance vector over time. Can use diff to generate the
%velocity vector. 
smoothVel = diff(smoothDist); %This ends up still having a lot of noise. Smooth out by 50 ms window
smoothVel = smooth(smoothVel,51);
smoothVel = smoothVel * 1000; %this is to produce cm/s

%reshape and make distance and time vectors into two column arrays
smoothDistTime = reshape(smoothDistTime,[],1);
smoothDist = reshape(smoothDist,[],1);
smoothVel = reshape(smoothVel,[],1);

disp('Finished Interpolation And Remapping of Velocity/Distance')



disp('Finished Finding Locomotion Starts')


funcOut=struct;
funcOut.Distance = [smoothDistTime,smoothDist];
funcOut.Velocity = [smoothVel];
% funcOut.BinaryLocomotion = locoBinary;
% funcOut.SimpleLocomotion = locoSimp;
% x.BinaryLocoRev = revBinary;
% funcOut.LocoStarts = smoothDistTime(locoStarts);
% funcOut.LocoEnds = smoothDistTime(locoEnds);


end


