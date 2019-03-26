%Test code to extract width by using a peak finding system followed by
%looking for X% of peak cutoffs to either side. Will search from peak out
%and from sides in. Will report two values for frequencies greater than and
%less than peak frequency respectively (one from outward search from peak,
%one from inward search from outside). 



function [widthVals,maxPos,maxVal,cutVal] = functionHeightBasedTuningWidth(tuningCurve,percentageCutoff,smoother);

%for testing
% smoother = 3;
% percentageCutoff = 0.10;% in percent


newCurve = smooth(tuningCurve,smoother);
%interpolate to 0.1 of the original sampling
newFreqVect = [1:0.1:length(newCurve)];
newCurve = interp1([1:1:length(newCurve)],newCurve,[1:0.1:length(newCurve)]);
% newCurve = interp1([1:1:length(newCurve)],newCurve,[1:0.1:length(newCurve)],'spline');

%now find peak
[maxVal,ind] = max(newCurve);

%need to deal with trash here. 
if maxVal <= 0
    widthVals = NaN(4,1);
    maxPos = NaN;
    maxVal = NaN;
    cutVal = NaN;
    disp('No Peak, All Values 0 or Negative')
else
    disp('Found Peak!')
    maxPos = newFreqVect(ind);
    cutVal = maxVal*percentageCutoff;
    %now focus below the peak value. 
    findLowOut = find(newCurve(1:ind) <= cutVal,1,'last');
    findLowIn = find(newCurve(1:ind) >= cutVal,1,'first');

    %now focus before peak value
    findHiOut = find(newCurve(ind:end) <= cutVal,1,'first');
    findHiIn = find(newCurve(ind:end) >= cutVal,1,'last');

    widthVals = NaN(4,1);
    if findLowOut
        widthVals(1) = newFreqVect(ind) - newFreqVect(findLowOut);
    end
    if findLowIn
        widthVals(2) = newFreqVect(ind) - newFreqVect(findLowIn);
    end
    if findHiOut
        widthVals(3) = newFreqVect(findHiOut+ind-1) - newFreqVect(ind);
    end
    if findHiIn
        widthVals(4) = newFreqVect(findHiIn+ind-1) - newFreqVect(ind);
    end
end




