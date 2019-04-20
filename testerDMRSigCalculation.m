%this is meant to be master code to run analysis code!
%Establishes the home folder as the base
masterFolder = pwd;
masterDir = dir;
masterDir = {masterDir.name};
%this code here is to remove folders reported in dir like "." or ".."
masterIndex = strfind(masterDir,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));


%masterFolders is a listing of all folders that I want to target for batch
%analysis
masterFolders = masterDir(masterIndex);

[findString] = functionCellStringFind(masterFolders,'.pdf');
masterFolders(findString) = [];



numFolders = length(masterFolders);
for masterCount = 1:numFolders
    disp('New File')
    %generates text that signifies path to the target folder
    targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
    cd(targetPath)

    fileName = strcat(masterFolders{masterCount},'DMRData.mat');
    load(fileName);
    
    numUnits = size(spikeArray,1);
    %now lets test significance of the responses. Since this requires many
    %loops to get things done, I think the best idea is to actually downsample
    %down to a 5ms time step, so that I can still get through things, show
    %significance 

    %the way in which I will do that is by smoothing by 5, then downsampling by
    %5. 
    dsLags = 20;
    dsStimulus = [];
    for i = 1:length(faxis)
        dsStimulus(i,:) = downsample(smooth(stimulus(i,:),5),5)*5;
    end
    
    dsSpikes = [];
    for i = 1:numUnits
        dsSpikes(i,:) = downsample(smooth(spikeArray(i,:),5),5)*5;
    end

    dsDMRtimes = downsample(spikeHistVect,5); %linspace(spikeHistVect(1),spikeHistVect(end),length(dsStimulus));

%     dsDMRtimes = interp1(ttlOnsetTime,dmrDIO,tempVect);
% 
%     dsDmrStep = mean(mode(diff(dsDMRtimes)));
%     % newWin = round(s.Parameters.STRFWin/dsDmrStep);
% 
%     dsSpikeHistVect = dsDMRtimes+dsDmrStep/2;

    %Now lets generate split STAs to calculate correlation. 
    randVals = randperm(length(dsSpikes),floor(length(dsSpikes)/2));
    splitSpikes1 = zeros(size(dsSpikes));
    splitSpikes1(:,randVals) = dsSpikes(:,randVals);
    randVals2 = [1:1:length(dsSpikes)];
    randVals2(randVals) = [];
    splitSpikes2 = zeros(size(dsSpikes));
    splitSpikes2(:,randVals2) = dsSpikes(:,randVals2);

    [sta1, stabigmat, spkcountvec] = quick_calc_sta(dsStimulus, splitSpikes1, dsLags);
    [sta2, stabigmat, spkcountvec] = quick_calc_sta(dsStimulus, splitSpikes2, dsLags);

    %now lets generate correlation coefficients. 
    for i = 1:numUnits
        tester = corrcoef(sta1(i,:),sta2(i,:));
        realCorrStore(i) = tester(2);
    end

    %now lets try and do the rotating shift of the spike train. 
    numShifts = 1000;
    shiftVal = floor(length(dsSpikes)/numShifts);

    corrCoeffStore= zeros(numShifts,numUnits);
    for i = 1:numShifts
        disp(strcat('Shift Number-',num2str(i)))
        %shuffle the spikes
        altArray =zeros(size(dsSpikes));
        altArray(:,1:length(dsSpikes) - shiftVal*i) = dsSpikes(:,shiftVal*i+1:end);
        altArray(:,length(dsSpikes) - shiftVal*i+1:length(dsSpikes)) = dsSpikes(:,1:shiftVal*i);
        %split in half? 
        tester = randperm((length(dsSpikes)),floor(length(dsSpikes)/2));
        altHalf1 = zeros(numUnits,length(dsSpikes));
        altHalf1(:,tester) = altArray(:,tester);
        revTest = [1:1:length(dsSpikes)];
        revTest(tester) = [];
        altHalf2 = zeros(numUnits,length(dsSpikes));
        altHalf2(:,revTest) = altArray(:,revTest);
        %now generate new STA
        [staAlt1, stabigmat, spkcountvec] = quick_calc_sta(dsStimulus, altHalf1, dsLags);
        [staAlt2, stabigmat, spkcountvec] = quick_calc_sta(dsStimulus, altHalf2, dsLags);
        %now we need to do unit by unit correlation coefficients
        for j = 1:numUnits
            tester = corrcoef(staAlt1(j,:),staAlt2(j,:));
            corrCoefStore(j,i) = tester(2);
        end
    end
    
    save(fileName,'faxis','stimulus','spikeHistVect','spikeArray','sta','realCorrStore','corrCoefStore')
    
    cd(masterFolder)
end



