%This function is meant to find the start and end of a response in order to
%determine the "width" of a response. The goal here will be to use a peak
%based detection method. Essentially, function will find the peak window
%within the proscribed region, determine what 15% of the difference of the
%peak and mean values are, and then use this to find the start and end. I
%will implement it as a crawling program so that it will detect dips in the
%middle of a response. 


function [widthOut] = functionResponseWidth(sigOut);


%variables I want to control
smoothFactor = 3;%how many bins to smooth by. 
percentCutoff = 0.15; %what percent of difference of peak vs baseline I want to use.
interpStep = 0.1; %factor of interpolation for signal

%pull out significance binary values

respSig = sigOut.Histogram(:,2);
bigSig = sigOut.Histogram(:,4);

timeVector = sigOut.Centers;
respHist = sigOut.Histogram(:,1);
smoothResp = smooth(respHist,smoothFactor);

timeStep = mean(diff(timeVector));
interpTime = [timeVector(1):timeStep*interpStep:timeVector(end)];
interpSmoothResp = interp1(timeVector,smoothResp,interpTime);

baselineHist = sigOut.BaselineHist;
baselineMean = mean(baselineHist);
%determine if response is single or multiple phase response
checkInd = 2;
whileInd = 0;
respNum = 0;
tempCount = 0;
widthStore = [];
startEndStore= [];
maxStore = [];
while whileInd == 0
%     tempCount
    if checkInd == length(respSig)
        disp('Reached End of Histogram')
        if tempCount > 0
            if tempCount <=2
                disp('Too small, not counting')
                tempCount = 0;
                checkInd = checkInd + 1;
                respNum = respNum - 1;
            else
                disp('large enough to count, check for big sig')
                endVal = checkInd; %now we have start and end values. 
                bigSigTest = mean(bigSig(startVal:endVal));
                if bigSigTest == 0
                    disp('No Big Sig, not Counting!')
                    tempCount = 0;
                    checkInd = checkInd + 1;
                    respNum = respNum - 1;
                else
                    tempCount = 0; %reset other values
                    checkInd = checkInd + 1; %advance checkInd
                    %now use start and end values to find peak
                    findStart = find(interpTime >= timeVector(startVal),1,'first');
                    findEnd = find(interpTime >= timeVector(endVal),1,'first');
                    [locMax maxInd] = max(interpSmoothResp(findStart:findEnd));
                    maxInd = maxInd + findStart - 1;
                    %find difference of baseline with peak
                    peakDiff = locMax - baselineMean;
                    %find X% of this difference, add to baseline mean
                    widthCutoff = baselineMean + peakDiff * percentCutoff;
                    %now find the points that demarcate the beginning and end of
                    %the response

                    %find width values
%                     maxInd
%                     length(interpSmoothResp)
                    respStart = find(interpSmoothResp(1:maxInd) <= widthCutoff,1,'last');
                    
                    disp('Found Start')
                    respEnd = find(interpSmoothResp(maxInd:length(interpSmoothResp)) <= widthCutoff,1,'first');
                    if length(respEnd) == 0
                        respEnd = length(interpSmoothResp) - maxInd;
                    end
                    disp('Found End')
                    widthStore(respNum) = interpTime(respEnd + maxInd) - interpTime(respStart);
                    startEndStore(:,respNum) = [interpTime(respStart),interpTime(respEnd + maxInd)];
                    maxStore(respNum) = interpTime(maxInd);
                    disp(strcat('Width Measured as:',num2str(widthStore(respNum))))
                end
            end
            break
        else
            break
        end
        break
    end
    %check value of currently indexed value
    prevVal = respSig(checkInd - 1);
    currVal = respSig(checkInd);
    if checkInd == 2 & prevVal == 1
        disp('ERROR: FIRST VALUE SIGNIFICANT, STARTING COUNT AS RESPONSE')
        startVal = 1;
        respNum = respNum + 1;
        tempCount = tempCount + 1;
        checkInd = checkInd + 1;
    elseif currVal - prevVal == 1
        disp('Detected Positive Shift, Starting Count!')
        startVal = checkInd;
        respNum = respNum + 1;
        tempCount = tempCount + 1;
        checkInd = checkInd + 1;
    elseif currVal - prevVal == -1
        disp('Detected Negative Shift, Ending Count!')
        %check the current count, if 1-2, delete
        if tempCount <=2
            disp('Too small, not counting')
            tempCount = 0;
            checkInd = checkInd + 1;
            respNum = respNum - 1;
        else
            disp('large enough to count, check for big sig')
            endVal = checkInd; %now we have start and end values. 
            bigSigTest = mean(bigSig(startVal:endVal));
            if bigSigTest == 0
                disp('No Big Sig, not Counting!')
                tempCount = 0;
                checkInd = checkInd + 1;
                respNum = respNum - 1;
            else
                disp('BigSig Checks out!')
%                 checkInd
                tempCount = 0; %reset other values
                checkInd = checkInd + 1; %advance checkInd
                %now use start and end values to find peak
                findStart = find(interpTime >= timeVector(startVal),1,'first');
                findEnd = find(interpTime >= timeVector(endVal),1,'first');
                [locMax maxInd] = max(interpSmoothResp(findStart:findEnd));
                maxInd = maxInd + findStart - 1;
                %find difference of baseline with peak
                peakDiff = locMax - baselineMean;
                %find X% of this difference, add to baseline mean
                widthCutoff = baselineMean + peakDiff * percentCutoff;
                %now find the points that demarcate the beginning and end of
                %the response

                %find width values
                respStart = find(interpSmoothResp(1:maxInd) <= widthCutoff,1,'last');
                
                
%                 figure
%                 plot(interpSmoothResp)
                if isempty(respStart) & startVal == 1
                    
                    respStart = 1;
                elseif isempty(respStart) & startVal ~= 1
                    respStart = find(interpTime >= timeVector(startVal),1,'first');
%                     startVal
%                     widthCutoff
%                     maxInd
%                     find(interpSmoothResp(1:maxInd) <= widthCutoff)
%                     error('FAILURE TO FIND START OF RESPONSE')
                end
%                 disp('Found Start')

                respEnd = find(interpSmoothResp(maxInd:length(interpSmoothResp)) <= widthCutoff,1,'first');
                if isempty(respEnd)
                    respEnd = length(interpSmoothResp) - maxInd - 1;
                end

%                 try
%                 disp('Found End')
                if (respEnd + maxInd) > length(interpSmoothResp)
                    widthStore(respNum) = interpTime(end) - interpTime(respStart);
                    startEndStore(:,respNum) = [interpTime(respStart),interpTime(end)];
                else
                    widthStore(respNum) = interpTime(respEnd + maxInd) - interpTime(respStart);
                    startEndStore(:,respNum) = [interpTime(respStart),interpTime(respEnd + maxInd)];
                end
                maxStore(respNum) = interpTime(maxInd);
%                 catch
%                     figure
%                     plot(interpSmoothResp)
%                 end
                
                
                disp(strcat('Width Measured as:',num2str(widthStore(respNum))))
            end
        end
        
    elseif currVal - prevVal == 0 & currVal == 1
        %this is simply continuation. 
        tempCount = tempCount + 1;
        checkInd = checkInd + 1;
%         disp('Continued Significance')
    elseif currVal - prevVal == 0 & currVal == 0
        checkInd = checkInd + 1;
%         disp('Continued Non-Significance')
    end
end


widthOut.Widths = widthStore;
widthOut.StartsEnds = startEndStore;
widthOut.PeakInds = maxStore;
disp('End Width Analysis')


end