%This is code that will read and link MBED inputs, photometry data, and
%noldus tracking information!
















%% First, lets pull the MBED stuff


[portStates] = maxTrialVariableNoTask(fname);

%we want to look at instates. Port 1 is TDT for photometry, port 2 is
%NOLDUS

inputPhot = [portStates.tStamps',portStates.inStates(:,1)];
inputNold = [portStates.tStamps',portStates.inStates(:,2)];

%now eliminate duplicates in photometry inputs
whileTrig = 1;
whileInd = 2;

while whileTrig == 1;
    if whileInd > length(inputPhot)
        break
    end
    prevVal = inputPhot(whileInd-1,2);
    currVal = inputPhot(whileInd,2);
    if prevVal == currVal
%         disp('Duplicate Detected!')
        inputPhot(whileInd,:) = [];
    else
        whileInd = whileInd + 1;
    end
end

%do the same for noldus inputs
whileTrig = 1;
whileInd = 2;

while whileTrig == 1;
    if whileInd > length(inputNold)
        break
    end
    prevVal = inputNold(whileInd-1,2);
    currVal = inputNold(whileInd,2);
    if prevVal == currVal
%         disp('Duplicate Detected!')
        inputNold(whileInd,:) = [];
    else
        whileInd = whileInd + 1;
    end
end

%pull just the onset times for photometry
onsetPhot = inputPhot(inputPhot(:,2) == 1,1);
onsetPhotDiff = diff(onsetPhot);

%clean up the signal
bigFinder = find(onsetPhotDiff>1050);
if length(bigFinder) == 1
    if bigFinder < 100
        onsetPhot(1:bigFinder) = [];
    elseif bigFinder > length(onsetPhot)-100
        onsetPhot(bigFinder:end) = [];
    end
elseif length(bigFinder)>1
    error('MORE THAN ONE FAILED PHOTOMETRY PULSE')
end

onsetPhotDiff = diff(onsetPhot);

%% Now lets pull the photometry inputs

%pull photometry trace
tracePhoto = data.streams.x70G.data;
%pull isobestic trace
traceIso = data.streams.x05G.data;
%pull timestamps for fluorescence
traceTiming = [0:1/data.streams.x70G.fs:(1/data.streams.x70G.fs)*(length(data.streams.x70G.data)-1)];
%pull jittered signal
traceJitt = data.epocs.PtE1.onset;
traceJittDiff = diff(traceJitt);

%check alignment
[xcf,lags,bounds]  = crosscorr(onsetPhotDiff,traceJittDiff);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);

if xcLag ~= 0
    error('Photometry Not Aligned')
elseif length(onsetPhot) ~= length(traceJitt)
    error('Mismatch in Number of Photometry Pulses')
end

traceDF = (tracePhoto - traceIso)./ traceIso;

































