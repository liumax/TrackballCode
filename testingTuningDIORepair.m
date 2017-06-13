%This is code to intelligently try and repair issues with ttls coming in
%from tuning

%This is stuff that is already set up
dioTimes = s.SoundData.ToneTimes;
dioTimeDiff = diff(dioTimes);

%pull this from s structure
predTimes = s.SoundData.Delays;


if length(predTimes) == length(dioTimes); %normal case
    disp('Times Matched Up with DIO Inputs')
elseif length(predTimes) ~= length(dioTimes); %case with fewer DIO times than predicted times
    disp(strcat('Mismatched DIO! Predicted:',num2str(length(predTimes)),',Actual:',num2str(length(dioTimes))))
    
    %now we trigger a while loop to look through and pull out differences
    %between predicted delays and actual delays. It looks like there is a
    %fairly consistent difference (about 0.17 seconds), but I dont think we
    %can rely on that. Instead, we will pick anything less than the minimum
    %predicted value. This will flag any things that are fucked and see if
    %it can remove them.
    
    %find minimum time between predicted ITIs. Any failure should exceed
    %this number in terms of the difference of the timing between predicted
    %and actual
    minPred = min(predTimes);
    lengthCheck = length(dioTimes) - length(predTimes);
    s.SoundData.BackupTrialMatrix = s.SoundData.TrialMatrix;
%     
%     %now do cross correlogram to determine whether the things are aligned
%     %properly. If not, then need to adjust accordingly
%     [xcf,lags,bounds]  = crosscorr(dioTimeDiff,predTimes);
%     %here, negative values mean that predicted samples are shifted left.
%     %This means that the first predicted sample is actually the second DIO
%     %sample.
%     [xcMax maxInd] = max(xcf);
%     xcLag = lags(maxInd);
%     
%     if xcLag == 0
%         disp('Samples are appropriately aligned')
%     elseif xcLag < 0
%         disp('Shifting ')
%     end
%     elseif xcLag < 0
%         disp('Predicted Samples Dont Capture All DIO. Removing Excess DIO')
%         dioTimes(1:xcLag) = [];
%         dioTimeDiff = diff(dioTimes);
%         disp('Shifted DIO times!')
%         s.SoundData.FlaggedTimes = [1:xcLag];
%         s.SoundData.FlagType = 'DIO Removal';
%     elseif xcLag > 0
%         disp('DIO is missing the correct ')
%     else
%     end
    repArray  = zeros(numFreqs,numDBs); %make array of zeros!
    alertArray = zeros(numFreqs,numDBs);
    %assume that the first TTL is fine
    freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(1,2));
    dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(1,3));
    repArray(freqFind,dbFind) = repArray(freqFind,dbFind) + 1;
    whileTrig = 1;
    whileCounter = 1;
    alertCounter = 1;
    while whileTrig == 1
        %first see if I've exceeded the bounds
        if whileCounter == length(predTimes) %if hit bound of predicted times]
            disp('Went through all predicted times')
            break
            
        elseif whileCounter == length(dioTimeDiff) %if hit bound of received times
            disp('Went through all received times')
            break
            
        else
            %check the difference between the current values. If exceeds
            %threshold, flag. Otherwise, need to record the tuning
            %information so I get accurate representation of repetitions
            diffCheck = dioTimeDiff(whileCounter) - predTimes(whileCounter);
            if diffCheck < minPred & diffCheck > 0
                %now lets fill in rep numbers
                freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(whileCounter + 1,2));
                dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(whileCounter + 1,3));
                repArray(freqFind,dbFind) = repArray(freqFind,dbFind) + 1;
                whileCounter = whileCounter + 1;
            elseif diffCheck >= minPred
                disp('Error Found! Insufficient DIO')
                freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(whileCounter + 1,2));
                dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(whileCounter + 1,3));
                %store information about the failed DIO
                alertArray(freqFind,dbFind) = repArray(freqFind,dbFind) + 1;
                alertDesig(alertCounter) = whileCounter;
                s.SoundData.TrialMatrix(whileCounter,:) = [];
                
                %dont need to update whileCounter, since this will cycle
                %through again.
            elseif diffCheck < 0
                disp('Error Found! Excess DIO')
                freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(whileCounter + 1,2));
                dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(whileCounter + 1,3));
                %delete the excess DIO
                dioTimes(whileCounter + 1) = [];
                dioTimeDiff = diff(dioTimes);
            else
                error('Code Failure In DIO Repair')
            end
        end
    end
end



for i = 1:10
    dioSample = dioTimeDiff(i:49+i);
    predSample = predTimes(1:50);
    testVal(i,1) = mean(dioSample - predSample);
    testVal(i,2) = std(dioSample-predSample);
end

