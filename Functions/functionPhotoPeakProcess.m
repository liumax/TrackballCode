
%This code is meant to process raw photometry data, smooth, pull peaks, and
%also remove decay. 

function [t_ds,newSmoothDS,targetPeaks] = functionPhotoPeakProcess(t,corrData,peakThresh);


%so first things first, we want to take the signal and remove the decay
%aspect. Fit a double exponential seems to work using the fit function. 

%fit to double exponential
f = fit(double(t)',corrData','exp2');
%plot fit
% figure
% plot(f,t,corrData);

%generate the fitted function, subtract
expFit = f.a*(exp(f.b*t)) + f.c*(exp(f.d*t));
newData = (corrData - expFit)./expFit;

%now we want to smooth. Since we are using 5ms bins, 21 bin smoothing should
%provide 50 ms smoothing in each direction. 

newSmooth = smooth(newData,21);

%since this runs at 200 hz (5ms bins), I can instead go down to 50 Hz,
%which is entirely reasonable.
t_ds = downsample(t,4);
newSmoothDS = downsample(newSmooth,4);


%now we need to pick out peaks. 

dSmoothDS = diff(newSmoothDS);
shifter = zeros(100,3);
shiftInd = 1;
pSign = 1;


for crawlInd = 1:length(dSmoothDS)
    %for the first point.
    if crawlInd == 1 & dSmoothDS(crawlInd) > 0
        pSign = 1;
    elseif crawlInd == 1 & dSmoothDS(crawlInd) < 0
        pSign = -1;
    elseif crawlInd == 1 & dSmoothDS(crawlInd) == 0
        whileTrig = 0;
        whilePlus = 1;
        while whileTrig == 0
            if dSmoothDS(crawlInd + whilePlus) >0
                shifter(shiftInd,1) = crawlInd;
                shifter(shiftInd,2) = 0;
                shifter(shiftInd,3) = 1;
                shiftInd = shiftInd + 1;
                whileTrig = 1;
            elseif dSmoothDS(crawlInd + whilePlus) <0
                shifter(shiftInd,1) = crawlInd;
                shifter(shiftInd,2) = 0;
                shifter(shiftInd,3) = -1;
                shiftInd = shiftInd + 1;
                whileTrig = 1;
            elseif dSmoothDS(crawlInd + whilePlus) == 0
                whilePlus = whilePlus + 1;
            end
        end
        %now for later data points
        %in the case of positive value
    elseif crawlInd >1 & dSmoothDS(crawlInd) > 0
        if pSign == -1
%             disp('Change Detected')
            shifter(shiftInd,1) = crawlInd;
            shifter(shiftInd,2) = pSign;
            shifter(shiftInd,3) = 1;
            shiftInd = shiftInd + 1;
            pSign = 1;
        elseif pSign == 0
%             disp('Change Detected')
            shifter(shiftInd,1) = crawlInd;
            shifter(shiftInd,2) = pSign;
            shifter(shiftInd,3) = 1;
            shiftInd = shiftInd + 1;
            pSign = 1;
        end
        %in the case of negative value
    elseif crawlInd > 1 & dSmoothDS(crawlInd) < 0
        if pSign == 1
%             disp('Change Detected')
            shifter(shiftInd,1) = crawlInd;
            shifter(shiftInd,2) = pSign;
            shifter(shiftInd,3) = -1;
            shiftInd = shiftInd + 1;
            pSign = -1;
        elseif pSign == 0
%             disp('Change Detected')
            shifter(shiftInd,1) = crawlInd;
            shifter(shiftInd,2) = pSign;
            shifter(shiftInd,3) = -1;
            shiftInd = shiftInd + 1;
            pSign = -1;
        end
    elseif crawlInd > 1 & dSmoothDS(crawlInd) == 0
        pSign = 0;
    end
end

%now go through this stuff

%classify peaks
for crawlInd = 1:length(shifter)
    if shifter(crawlInd,2) == -1 & shifter(crawlInd,3) == 1
        shifter(crawlInd,4) = 1; %This is the dip
    elseif shifter(crawlInd,2) == 1 & shifter(crawlInd,3) == -1
        shifter(crawlInd,4) = 2; %This is the peak
    elseif shifter(crawlInd,2) == 0 & shifter(crawlInd,3) == 1
        shifter(crawlInd,4) = 0; % positive inflection
    elseif shifter(crawlInd,2) == 0 & shifter(crawlInd,3) == -1
        shifter(crawlInd,4) = 0; %negative inflection
    elseif shifter(crawlInd,2) == 1 & shifter(crawlInd,3) == 0
        shifter(crawlInd,4) = 0; %kind of a trough?
    elseif shifter(crawlInd,2) == -1 & shifter(crawlInd,3) == 0
        shifter(crawlInd,4) = 0; %kind of a peak?
    end
end

%basically, I want to eliminate any peaks that occur with no negative peak
%before

whileTrig = 0;

while whileTrig == 0;
    %check if first shifter value is 1. 
    if shifter(1,4) == 2;
        shifter(1,:) = [];
        disp('Removing Peak with No Trough')
    elseif shifter(1,4) ~= 2;
        whileTrig = 1;
    end
end

%now go through and eliminate any things where there is no state change (no
%change in values

whileTrig = 0;
whileCounter = 2;
prevValue = shifter(1,4);
while whileTrig == 0;
    if whileCounter == length(shifter)
        whileTrig = 1;
    end
    currValue = shifter(whileCounter,4);
    if currValue == prevValue
        shifter(whileCounter,:) = [];
        disp('Cutting Duplicate')
    else
        prevValue = currValue;
        whileCounter = whileCounter + 1;
    end
end

%now remove the last negative peak, so there is only a positive peak at the
%end
if shifter(end,4) ~= 2
    shifter(end,:) = [];
    disp('Trimming Ends')
end


%now get peak and trough values.
for crawlInd = 1:length(shifter)
    shifter(crawlInd,5) = newSmoothDS(shifter(crawlInd,1));
end

peakInds = find(shifter(:,4) == 2);

peakVals(:,1) = shifter(peakInds,5) - shifter(peakInds-1,5);
peakVals(:,2) = shifter(peakInds,5);
peakVals(:,3) = shifter(peakInds-1,5);

%now, sort peaks by size. 
[Y I] = sort(peakVals(:,1));

%so now we want to institute some kind of size and timing exclusion? 

%actually, we arent going to do that. Just do size exclusion. Set limit to
%0.01

% peakThresh = 0.01;
targetPeaks = peakVals(peakVals(:,1) > peakThresh,:);
% 
% figure
% plot(t_ds,newSmoothDS,'k')
% hold on
% plot(t_ds(targetPeaks(:,2)),newSmoothDS(targetPeaks(:,2)),'r*')
% plot(t_ds(targetPeaks(:,3)),newSmoothDS(targetPeaks(:,3)),'g*')
% xlim([t_ds(1) t_ds(end)])
% ylabel('dF/F')
% xlabel('Time (s)')
% title('Peak Selection')


end























