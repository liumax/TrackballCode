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
rewList = zeros(10,2);
rewCounter = 1;

for i=1:length(fileNames)
    fileName = fileNames{i};
    %extract data
    [s] = functionPhotoTwoTonePavExtract(fileName,[-3,5],0.01,0);
    %extract number of trials that things have been updated by
    trialList(i,1) = trialCounter;
    trialList(i,2) = length(s.Trials);
    disp('updated trial list')
    
    %extract number of rewarded trials
    rewList(i,1) = rewCounter;
    rewList(i,2) = length(find(s.Trials(:,2)>0));
    disp('updated reward list')
    %now lets put data into bigger arrays
    
    tonePhot(:,trialCounter:trialCounter + trialList(i,2)-1) = s.Photometry.ToneRaster;
    rewPhot(:,rewCounter:rewCounter + rewList(i,2)-1) = s.Photometry.RewardRaster;
%     
%     toneLick() = ;
%     rewLick() = ;
%     
%     toneVel = ;
%     rewVel = ;
%     
    trialCounter = trialCounter + trialList(i,2);
    rewCounter = rewCounter + rewList(i,2);
end





