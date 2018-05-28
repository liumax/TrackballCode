%this code is meant to take a simplified locomotion trace from the rotary
%encoder (+1 for positive movement, -2 for negative movement), and process
%it to a.) clean up ends, and b.) merge periods together. This is agnostic
%to negative or positive, and should be capable of working with either. 

% inputs: 
% 
% locoSimp (1 column vector)
% threshLim (limit for merging loco bouts)
% interpStep (to normalize threshold limit)
% values to pick up:val1 val2

function [locoOut] = locoBinFix(locoSimp,threshLim,interpStep,val1,val2);

locoStarts = find(diff(locoSimp) == val1 | diff(locoSimp) == val2);
locoEnds = find(diff(locoSimp) == -val1| diff(locoSimp) == -val2)+1;


disp('Finding Locomotion Starts')


%figure out if a start or an end comes first
[minVal minInd] = min([locoStarts(1),locoEnds(1)]);

%if locomotion starts are first, this is fine, and as expected. 
%if locomotion ends come first, we need to create a false start at the
%beginning
if minInd == 2
    locoStarts(2:end+1) = locoStarts(1:end);
    locoStarts(1) = 1;
    disp('Animal began while locomoting, Added First Value to LocoStarts')
else
    disp('Animal Starts Stationary')
end

%now lets check to see if we have equal numbers of starts and stops. 
if length(locoStarts) == length(locoEnds)
    disp('Same number of locomotion starts and ends')
else
    disp('Different number of locomotion starts and ends...examining...')
    %only way this should happen should be if animal is still
    %locomoting at the end of the session. Check if last value is
    %positive
    if locoSimp(end) > 0
        disp('Animal Locomoting through end of recording, adding fake stop at end')
        locoEnds(end+1) = length(locoSimp);
    else
        error('INCORRECT NUMBER OF STARTS AND STOPS FOR FORWARD LOCOMOTION')
    end
end

%now that inputs and outputs should have been fixed, lets try and actually
%fix spacing. 
loopTrig = 0;
% threshLim = 1;
loopInd = 1;
newLocoSimp = locoSimp;
    
%convert threshold to the index of locomotion time steps
threshVal = threshLim/interpStep;
while loopTrig == 0
    %now we go through and advance single locomotion bout at a time. Idea
    %is that if distance from nearest end to next start doesnt exceed
    %threshold limit, then merge the two. 
    if loopInd == length(locoStarts)
        disp('Finished Start/End Merging')
        break
    end
    currEnd = locoEnds(loopInd);
    nextStart = locoStarts(loopInd+1);
    if nextStart - currEnd < threshVal
%         disp('Small Timing Detected')
        %get velocity value in between. If == 0, then proceed. 
        if mean(locoSimp(currEnd:nextStart)) == 0;
            locoStarts(loopInd + 1) = [];
            locoEnds(loopInd) = [];
%             disp(strcat('Merged',num2str(loopInd),'&',num2str(loopInd+1)))
            newLocoSimp(currEnd:nextStart) = val1;
        else
%             disp('Found Non-Zero Velocity In between, cancel merge')
            loopInd = loopInd + 1;
        end
    else
        loopInd = loopInd + 1;
    end
end

locoOut.newSimp = newLocoSimp;
locoOut.starts = locoStarts;
locoOut.ends = locoEnds;















