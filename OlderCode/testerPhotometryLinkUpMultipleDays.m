%This is test code to link photometry data from multiple days, as well as
%trials, and process and display it together.

%load filenames as a cell array
fileNames = { 'ML170511ATaskNoRew170615','ML170511A_170619_TwoTonePavlov',...
    'ML170511A_170620_TwoTonePavlov','ML170511A_170621_TwoTonePavlov',...
    'ML170511A_170622_TwoTonePavlov','ML170511A_170623_TwoTonePavlov'};

%now we want a big for loop that goes through and pulls everything. 
mainStruct = struct;
trialList = zeros(10,2);
trialCounter = 1;
rewList = zeros(100,2);
rewCounter = 1;
lickStore = zeros(100,1);


for i=1:length(fileNames)
    fileName = fileNames{i};
    fileName = strcat(fileName,'Analysis');
    load(fileName);
    %extract data
    
    %pull number of trials, store.
    trialList(i,1) = trialCounter;
    trialList(i,2) = length(s.MBED.ToneDelivery);
    
    %extract identity of rewarded trials
    rewList(rewCounter:rewCounter + length(s.MBED.HiTrials) - 1,1) = s.MBED.HiTrials + trialCounter-1;
    rewList(rewCounter:rewCounter + length(s.MBED.HiTrials) - 1,2) = i;
    %extract traces
    tonePhot(:,trialCounter:trialCounter + trialList(i,2)-1) = s.PhotoRaster.ToneRaster;
    
        
    %extract licking times.
    for j = 1:length(s.MBED.HiTrials)
        findLicks = find(s.Licking.ToneRaster(:,2) == s.MBED.HiTrials(j) & s.Licking.ToneRaster(:,1) > 0 & s.Licking.ToneRaster(:,1) < 1);
        if length(findLicks)>0
            lickStore(rewCounter+j-1,1) = length(findLicks);
            
        end
        lickStore(rewCounter+j-1,2) = i;
    end
    
    trialCounter = trialCounter + trialList(i,2);
    rewCounter = rewCounter + length(s.MBED.HiTrials);

    
end



%% test code 170715

%divide rewarded and unrewarded trials into bins of 10 trials
tester = round(length(rewList)/10);

for i = 1:tester
    photStore(:,i) = mean(tonePhot(:,rewList(((i-1)*10)+1:(i*10))),2);
    %normalize based on first 2000 things
    preAv = mean(photStore(1:2000,i));
    photStore(:,i) = photStore(:,i) - preAv;
end

figure
imagesc(photStore')

%hmm try with individual trials

for i = 1:length(rewList)
    photStoreInd(:,i) = tonePhot(:,rewList(i));
    preAv = mean(photStoreInd(1:2000,i));
    photStoreInd(:,i) = photStoreInd(:,i) - preAv;
end

figure
imagesc(photStoreInd')

%Lets try this with unrewarded trials
allTrials = [1:size(tonePhot,2)];
noRew = setDiff(allTrials,rewList);

tester = round(length(noRew)/10);

for i = 1:tester
    photStore(:,i) = mean(tonePhot(:,noRew(((i-1)*10)+1:(i*10))),2);
    %normalize based on first 2000 things
    preAv = mean(photStore(1:2000,i));
    photStore(:,i) = photStore(:,i) - preAv;
end

figure
imagesc(photStore')

%hmm try with individual trials
photStoreInd = [];
for i = 1:length(noRew)
    photStoreInd(:,i) = tonePhot(:,noRew(i));
    preAv = mean(photStoreInd(1:2000,i));
    photStoreInd(:,i) = photStoreInd(:,i) - preAv;
end

figure
imagesc(photStoreInd')

%lets look at transition from the first day to second day

photStoreInd = [];

for i = 90:110
    photStoreInd(:,i) = tonePhot(:,rewList(i));
    preAv = mean(photStoreInd(1:2000,i));
    photStoreInd(:,i) = photStoreInd(:,i) - preAv;
end

figure
imagesc(photStoreInd')


photStoreInd = [];

for i = 80:120
    photStoreInd(:,i) = tonePhot(:,noRew(i));
    preAv = mean(photStoreInd(1:2000,i));
    photStoreInd(:,i) = photStoreInd(:,i) - preAv;
end

figure
imagesc(photStoreInd')



%bumping each trace up so that I can see all of them

for i = 1:tester
    photStore(:,i) = mean(tonePhot(:,rewList(((i-1)*10)+1:(i*10))),2) - 0.02*(i-1);
end

