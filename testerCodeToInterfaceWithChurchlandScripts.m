%Lets make some code that is meant to convert a raster dataset of mine into
%the input for the fano factor code.

%fano factor code takes in a spike train block of data, which is an m x n
%matrix. m trials, n milliseconds, with all points being either ones or
%zeros. These are stored in a structured array in which there is one field
%(spikes), and there are rows based on different units

homeDir = '/Users/maxliu/testy/ChurchlandVariance/Variance_toolbox';

cd /Users/maxliu/Downloads/170414LocalAnalysis/DatStimControl
load('datStimControl64ChanFullData.mat')
fieldNames = fields(fullData);

spikeData = struct;

unitCounter = 1;

for i = 1:length(fieldNames);
    %open dataset
    testData = fullData.(fieldNames{i});
    numUnits = length(testData.DesignationName);
    rasterWindow = testData.Parameters.RasterWindow;
    numTrials = size(testData.(testData.DesignationName{1}).IndividualHistograms,2);
    %prep vector
    
    for j = 1:numUnits;
        fireVector = zeros(numTrials,(rasterWindow(2) - rasterWindow(1))*1000);
        rasterData = testData.(testData.DesignationName{j}).AllRasters;
        rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
        %remove points at which raster data == 0
        rasterData(rasterData(:,1) ==0,:) = [];
        for k = 1:numTrials
            spikeFind = find(rasterData(:,2) == k);
            if find(rasterData(:,2) == k) %if data exists for this time point
                fireVector(k,rasterData(spikeFind,1)) = 1;
            end
        end
        spikeData(unitCounter).spikes = logical(fireVector);
        unitCounter = unitCounter + 1;
    end
end


%set things for churchland code
times = 3800:25:4500; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%now lets try body stim
clear
cd /Users/maxliu/Downloads/170414LocalAnalysis/DatStimBody
load('170413PrunedBodyStimDataWith170413Data.mat')

fieldNames = fields(fullData);

spikeData = struct;

unitCounter = 1;

for i = 1:length(fieldNames);
    %open dataset
    testData = fullData.(fieldNames{i});
    numUnits = length(testData.DesignationName);
    rasterWindow = testData.Parameters.RasterWindow;
    numTrials = size(testData.(testData.DesignationName{1}).IndividualHistograms,2);
    %prep vector
    
    for j = 1:numUnits;
        fireVector = zeros(numTrials,(rasterWindow(2) - rasterWindow(1))*1000);
        rasterData = testData.(testData.DesignationName{j}).AllRasters;
        rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
        %remove points at which raster data == 0
        rasterData(rasterData(:,1) ==0,:) = [];
        for k = 1:numTrials
            spikeFind = find(rasterData(:,2) == k);
            if find(rasterData(:,2) == k) %if data exists for this time point
                fireVector(k,rasterData(spikeFind,1)) = 1;
            end
        end
        spikeData(unitCounter).spikes = logical(fireVector);
        unitCounter = unitCounter + 1;
    end
end


%set things for churchland code
times = 3800:25:5500; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 100;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%dont see much despite a bit of wiggling with settings
%% try with some sensory stimuli
clear
cd /Users/maxliu/Downloads
load('secondPairing.mat')


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tuning1.ToneDur;
numTrials = size(s.(s.DesignationName{1}).Tuning1Analysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tuning1Analysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 50:25:650; %this is centered around zero point
fanoParams.alignTime = 400; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);
%looks like i've got a suppression! Yay. 

%now look at the after

spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tuning1.ToneDur;
numTrials = size(s.(s.DesignationName{1}).Tuning2Analysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tuning2Analysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 50:25:650; %this is centered around zero point
fanoParams.alignTime = 400; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%this looks qualitatively similar. Maybe weaker? It looks like overall
%response has changed too

%Try looking at the pairing trials?

spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = [-1 2]; %here i know this is 1 second 
numTrials = size(s.(s.DesignationName{1}).PairingAnalysis.FirstSpikeTimes,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).PairingAnalysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 500:25:1500; %this is centered around zero point
fanoParams.alignTime = 1000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%see big reduction here! lets try and do things with first vs second

spikeData1 = struct;
spikeData2 = struct;
spikeData3 = struct;
spikeData4 = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = [-1 2]; %here i know this is 1 second 
numTrials = size(s.(s.DesignationName{1}).PairingAnalysis.FirstSpikeTimes,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).PairingAnalysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData1(unitCounter).spikes = logical(fireVector(1:50,:));
    spikeData2(unitCounter).spikes = logical(fireVector(50:100,:));
    spikeData3(unitCounter).spikes = logical(fireVector(101:150,:));
    spikeData4(unitCounter).spikes = logical(fireVector(151:200,:));
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 500:25:1500; %this is centered around zero point
fanoParams.alignTime = 1000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData1, times, fanoParams);
plotFano(Result);
Result = VarVsMean(spikeData2, times, fanoParams);
plotFano(Result);
Result = VarVsMean(spikeData3, times, fanoParams);
plotFano(Result);
Result = VarVsMean(spikeData4, times, fanoParams);
plotFano(Result);

%% now try looking at pairing data with the control/target

clear
load('170404DATPairingDataMinusDuplicate.mat')
s = fullData.x170223_ML170108A_R17_3000_TruePairingC4T8laser12mW;


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tuning1.ToneDur;
numTrials = size(s.(s.DesignationName{1}).Tuning1Analysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tuning1Analysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 50:25:650; %this is centered around zero point
fanoParams.alignTime = 400; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tuning1.ToneDur;
numTrials = size(s.(s.DesignationName{1}).Tuning2Analysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tuning2Analysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 50:25:650; %this is centered around zero point
fanoParams.alignTime = 400; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);
%odnt really see any drop in the fano. 

%% look at tone presentations before pairing

spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).Tones1ControlAnalysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tones1ControlAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 3500:100:4500; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).Tones1ControlAnalysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tones1ControlAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end

%set things for churchland code
times = 3500:100:4500; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%not seeing much!

%%  Try pairing
spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).PairingControlAnalysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).PairingControlAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).PairingTargetAnalysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).PairingTargetAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end

%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%% Try with the otehr dataset? 
clear
load('170404DATPairingDataMinusDuplicate.mat')
s = fullData.x170303_ML170108B_L17_3100_pairingT16C8db100;


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tuning1.ToneDur;
numTrials = size(s.(s.DesignationName{1}).Tuning1Analysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tuning1Analysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 50:25:650; %this is centered around zero point
fanoParams.alignTime = 400; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tuning1.ToneDur;
numTrials = size(s.(s.DesignationName{1}).Tuning2Analysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tuning2Analysis.AllRasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 50:25:650; %this is centered around zero point
fanoParams.alignTime = 400; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);
%In the second plot, see a complete and total drop out of the fano. Maybe a
%thing comes on before tune???? confusing, and code throws error : A is rank deficient to within machine precision. 

%% look at tone presentations before pairing

spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).Tones1ControlAnalysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tones1ControlAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 100;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).Tones1TargetAnalysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tones1TargetAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end

%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 100;
fanoParams.matchReps = 0; %had to add this to fix bugs


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%control vs target, seeing something of an increase with the target
%specifically??

%%  Try pairing
spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).PairingControlAnalysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).PairingControlAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;
fanoParams = rmfield(fanoParams,'matchReps'); %removing matchReps!


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones1.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).PairingTargetAnalysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).PairingTargetAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end

%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%In target, seeing a big increase following tone. 

%% Try post
spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones2.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).Tones2ControlAnalysis.IndivHistograms,2);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tones2ControlAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end


%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);


spikeData = struct;

unitCounter = 1;

numUnits = length(s.DesignationName);
rasterWindow = s.Parameters.rasterWindow*s.SoundData.Tones2.ToneDuration;
numTrials = size(s.(s.DesignationName{1}).Tones2TargetAnalysis.IndivHistograms,1);
%prep vector

for j = 1:numUnits;
    fireVector = zeros(numTrials,round((rasterWindow(2) - rasterWindow(1))*1000));
    rasterData = s.(s.DesignationName{j}).Tones2TargetAnalysis.Rasters;
    rasterData(:,1) = round((rasterData(:,1) - rasterWindow(1))*1000);
    %remove points at which raster data == 0
    rasterData(rasterData(:,1) ==0,:) = [];
    for k = 1:numTrials
        spikeFind = find(rasterData(:,2) == k);
        if find(rasterData(:,2) == k) %if data exists for this time point
            fireVector(k,rasterData(spikeFind,1)) = 1;
        end
    end
    spikeData(unitCounter).spikes = logical(fireVector);
    unitCounter = unitCounter + 1;
end

%set things for churchland code
times = 3000:100:6000; %this is centered around zero point
fanoParams.alignTime = 4000; % this time will become zero time 
fanoParams.boxWidth = 50;


Result = VarVsMean(spikeData, times, fanoParams);

plotFano(Result);

%hmmm...firing more or less seems to have gone away. but the response also
%seems to be gone?