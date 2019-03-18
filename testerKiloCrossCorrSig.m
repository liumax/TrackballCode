%This is a snippet of code meant to extract significant responses in the
%cross correlogram. 



for i = 1:numUnits
    for j = 1:numUnits
        %check to see if there is a jittered cross corr 
        if i > size(s.CrossCorrJitter,1) | j > size(s.CrossCorrJitter,2)
            disp('Reaching End of Jittered Corr')
            break
        end
        %now we need to determine if the cell is full
        if ~isempty(s.CrossCorrJitter{i,j})
            %pull the jittered values
            jitVals = s.CrossCorrJitter{i,j};
            %pull the actual cross corr
            trueVals = s.CrossCorrs{i,j};
            %pull percentiles
            Y = prctile(jitVals,[.5 99.5],2);
            %now lets make sure anything lower gets picked up
            sigVals = zeros(length(trueVals),1);
            sigVals(trueVals - Y(:,1)' < 0) = -1;
            %now do higher
            sigVals(trueVals - Y(:,2)' > 0) = 1;
            %now we want to clean this up. Remove all single change
            %stretches, only keep prolonged changes. 
            whileTrig = 0;
            creepCount = 1;
            minWidth = 2;
            while whileTrig == 0
                if creepCount > length(sigVals)
                    break
                end
                %check if sigVals == 0
                if sigVals(creepCount) == 0
                    creepCount = creepCount + 1;
                else
                    tmpVal = sigVals(creepCount);
                    storeCount = creepCount;
                    while whileTrig == 0
                        if tmpVal == sigVals(creepCount) %this is still part of a continuous grouping
                            creepCount = creepCount + 1;
                        elseif tmpVal ~= sigVals(creepCount) %no longer part of continuous block of numbers
                            %determine length
                            tmpLength = creepCount - storeCount;
                            if tmpLength > minWidth;
                                disp('Continuous Sig Detected')
                                break
                            else
                                disp('Insufficient Length, deleting...')
                                sigVals(storeCount:creepCount - 1) = 0;
                                break
                            end
                        end
                    end
                end
            end
            newStore{i,j} = sigVals;
            sigWarnPos(i,j) = length(find(sigVals == 1));
            sigWarnNeg(i,j) = length(find(sigVals == -1));
        end
    end
end