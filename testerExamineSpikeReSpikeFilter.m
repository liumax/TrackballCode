%this is code to validate the stuffs!Checking if I can pick spikes back out
%again based on timing. 
cd E:\170420CodeToTestReSortNewSpikeFilter\Comparison

%pull out the multiunit data
load('MUData.mat')
mudata = s;
s = [];

%pull up sorted dataset
load('170330_ML170321C_R05_2700_manualDatStimDatStimAnalysis.mat')

%lets try this in a singular case first

testData = s.nt2cluster1.SpikeTimes;
dataMU = mudata.nt2cluster1.SpikeTimes;
spikeLim = [-0.001 0.001];
zeroFind = zeros(length(testData),1);

counter = 1;
for i = 1:length(testData)
    testPoint = testData(i);
    spikeSub = dataMU - testPoint;
    %select points within certain limit
    spikeFinder = find(spikeSub > spikeLim(1) & spikeSub < spikeLim(2));
    spikeStore(counter: counter + length(spikeFinder) - 1,1) = spikeFinder;
    spikeStore(counter: counter + length(spikeFinder) - 1,2) = i;
    counter = counter + length(spikeFinder);
    spikeNum(i) = length(spikeFinder);
    if find(spikeSub == 0)
        zeroFind(i) = 1;
    end
end

%IT WORKS :DDDD
%try doing for another. 

testData = s.nt4cluster3.SpikeTimes;
dataMU = mudata.nt4cluster1.SpikeTimes;
spikeLim = [-0.001 0.001];
zeroFind = zeros(length(testData),1);

counter = 1;
for i = 1:length(testData)
    testPoint = testData(i);
    spikeSub = dataMU - testPoint;
    %select points within certain limit
    spikeFinder = find(spikeSub > spikeLim(1) & spikeSub < spikeLim(2));
    spikeStore(counter: counter + length(spikeFinder) - 1,1) = spikeFinder;
    spikeStore(counter: counter + length(spikeFinder) - 1,2) = i;
    counter = counter + length(spikeFinder);
    spikeNum(i) = length(spikeFinder);
    if find(spikeSub == 0)
        zeroFind(i) = 1;
    end
end
%works!


%Try for all units in the dataset
overCounter = 1;
%try for all
for j = 1:length(mudata.DesignationName)
    muData = mudata.(mudata.DesignationName{j}).SpikeTimes;
    %now I need to extract the names!
    muName = mudata.DesignationName{j};
    stringFind = strfind(muName,'cluster');
    trodeName = muName(1:stringFind);
    %now go to the clustered data and pull things out!
    IndexC = strfind(s.DesignationName, trodeName);
    Index = find(not(cellfun('isempty', IndexC)));
    for k = 1:length(Index)
        testData = s.(s.DesignationName{Index(k)}).SpikeTimes;
        zeroFind = zeros(length(testData),1);
        counter = 1;
        for i = 1:length(testData)
            testPoint = testData(i);
            spikeSub = dataMU - testPoint;
            %select points within certain limit
            spikeFinder = find(spikeSub > spikeLim(1) & spikeSub < spikeLim(2));
            spikeStore(counter: counter + length(spikeFinder) - 1,1) = spikeFinder;
            spikeStore(counter: counter + length(spikeFinder) - 1,2) = i;
            counter = counter + length(spikeFinder);
            spikeNum(i) = length(spikeFinder);
            if find(spikeSub == 0)
                zeroFind(i) = 1;
            end
        end
        spikeMatch(overCounter) = sum(zeroFind)/length(s.(s.DesignationName{Index(k)}).SpikeTimes);
        overCounter = overCounter + 1;
        cellMatch{overCounter,1} = muName;
        cellMatch{overCounter,2} = s.DesignationName{Index(k)};
    end
    
    %try and match the thing to the MU
    
    
end
%it turns out perfect matches are kind of hard...

%lets try this with narrower windows
overCounter = 1;
spikeLim = [-0.0005 0.0005];
%try for all
for j = 1:length(mudata.DesignationName)
    muData = mudata.(mudata.DesignationName{j}).SpikeTimes;
    %now I need to extract the names!
    muName = mudata.DesignationName{j};
    stringFind = strfind(muName,'cluster');
    trodeName = muName(1:stringFind);
    %now go to the clustered data and pull things out!
    IndexC = strfind(s.DesignationName, trodeName);
    Index = find(not(cellfun('isempty', IndexC)));
    for k = 1:length(Index)
        testData = s.(s.DesignationName{Index(k)}).SpikeTimes;
        zeroFind = zeros(length(testData),1);
        counter = 1;
        for i = 1:length(testData)
            testPoint = testData(i);
            spikeSub = dataMU - testPoint;
            %select points within certain limit
            spikeFinder = find(spikeSub > spikeLim(1) & spikeSub < spikeLim(2));
            spikeStore(counter: counter + length(spikeFinder) - 1,1) = spikeFinder;
            spikeStore(counter: counter + length(spikeFinder) - 1,2) = i;
            counter = counter + length(spikeFinder);
            spikeNum(i) = length(spikeFinder);
            if find(spikeSub == 0)
                zeroFind(i) = 1;
            end
        end
        spikeMatch(overCounter) = sum(zeroFind)/length(s.(s.DesignationName{Index(k)}).SpikeTimes);
        overallSpikeMatch(overCounter) = length(spikeStore)/length(s.(s.DesignationName{Index(k)}).SpikeTimes);
        overCounter = overCounter + 1;
        cellMatch{overCounter,1} = muName;
        cellMatch{overCounter,2} = s.DesignationName{Index(k)};
    end
    
    %try and match the thing to the MU
    
    
end

%This approach seems to have overproduced spikes for a number of cells, and
%underproduced for some others. Perhaps plotting histograms of these values
%would be helpful?
overCounter = 1;
spikeLim = [-0.0015 0.0015];
%try for all
for j = 1:length(mudata.DesignationName)
    muData = mudata.(mudata.DesignationName{j}).SpikeTimes;
    %now I need to extract the names!
    muName = mudata.DesignationName{j};
    stringFind = strfind(muName,'cluster');
    trodeName = muName(1:stringFind);
    %now go to the clustered data and pull things out!
    IndexC = strfind(s.DesignationName, trodeName);
    Index = find(not(cellfun('isempty', IndexC)));
    for k = 1:length(Index)
        testData = s.(s.DesignationName{Index(k)}).SpikeTimes;
        zeroFind = zeros(length(testData),1);
        counter = 1;
        for i = 1:length(testData)
            testPoint = testData(i);
            spikeSub = muData - testPoint;
            %select points within certain limit
            spikeFinder = find(spikeSub > spikeLim(1) & spikeSub < spikeLim(2));
            spikeStore(counter: counter + length(spikeFinder) - 1,1) = spikeSub(spikeFinder);
            spikeStore(counter: counter + length(spikeFinder) - 1,2) = i;
            counter = counter + length(spikeFinder);
            spikeNum(i) = length(spikeFinder);
            if find(spikeSub == 0)
                zeroFind(i) = 1;
            end
        end
        spikeMatch(overCounter) = sum(zeroFind)/length(s.(s.DesignationName{Index(k)}).SpikeTimes);
        overallSpikeMatch(overCounter) = length(spikeStore)/length(s.(s.DesignationName{Index(k)}).SpikeTimes);
        overCounter = overCounter + 1;
        cellMatch{overCounter,1} = muName;
        cellMatch{overCounter,2} = s.DesignationName{Index(k)};
        figure
        hist(spikeStore(:,1),[spikeLim(1):0.000001:spikeLim(2)])
        spikeStore = [];
    end
    
    %try and match the thing to the MU
    
    
end

%LOL TURNED OUT I HAD A BUGS
%Zero finding is good for this entire dataset it turns out. bwahahahahaha


%try with different dataset

clear
%pull out the multiunit data
load('MUActivity170412.mat')
mudata = s;
s = [];

%pull up sorted dataset
load('170412_ML170306C_R10_2500_datStimDatStimAnalysis.mat')


overCounter = 1;
spikeLim = [-0.0009 0.0009];
%try for all
for j = 1:length(mudata.DesignationName)
    muData = mudata.(mudata.DesignationName{j}).SpikeTimes;
    %now I need to extract the names!
    muName = mudata.DesignationName{j};
    stringFind = strfind(muName,'cluster');
    trodeName = muName(1:stringFind);
    %now go to the clustered data and pull things out!
    IndexC = strfind(s.DesignationName, trodeName);
    Index = find(not(cellfun('isempty', IndexC)));
    for k = 1:length(Index)
        testData = s.(s.DesignationName{Index(k)}).SpikeTimes;
        zeroFind = zeros(length(testData),1);
        counter = 1;
        for i = 1:length(testData)
            testPoint = testData(i);
            spikeSub = muData - testPoint;
            %select points within certain limit
            spikeFinder = find(spikeSub > spikeLim(1) & spikeSub < spikeLim(2));
            spikeStore(counter: counter + length(spikeFinder) - 1,1) = spikeSub(spikeFinder);
            spikeStore(counter: counter + length(spikeFinder) - 1,2) = i;
            counter = counter + length(spikeFinder);
            spikeNum(i) = length(spikeFinder);
            if find(spikeSub == 0)
                zeroFind(i) = 1;
            end
        end
        spikeMatch(overCounter) = sum(zeroFind)/length(s.(s.DesignationName{Index(k)}).SpikeTimes);
        overallSpikeMatch(overCounter) = length(spikeStore)/length(s.(s.DesignationName{Index(k)}).SpikeTimes);
        overCounter = overCounter + 1;
        cellMatch{overCounter,1} = muName;
        cellMatch{overCounter,2} = s.DesignationName{Index(k)};
        figure
        hist(spikeStore(:,1),[spikeLim(1):0.000001:spikeLim(2)])
        spikeStore = [];
    end
    
    %try and match the thing to the MU
    
    
end


%didnt work with 1.5 ms windows. Try 0.6 ms
%0.6 ms too tight. try 0.9






