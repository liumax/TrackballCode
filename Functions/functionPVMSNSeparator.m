
%This function takes in the average wave forms, selects the biggest one,
%and delivers a value (1 = PV 0 = MSN) based on the peak valley threshold,
%which I have temporarily set to 0.5 ms
function [funcOut,peakTrough] = functionPVMSNSeparator(waveForms,thresholdPeakValley);

%find size of waveforms
waveSize = size(waveForms);

%find maximum waveform by max peak size
maxWave = max(waveForms);
[maxVal maxInd] = max(maxWave);

%chose the big wave, interpolate to fine degree
chosenWave = waveForms(:,maxInd);
interpVect = [1:0.1:40];
interpWave = interp1(1:40,chosenWave,interpVect,'spline');

%now we need to find the peak. Find this starting at point 10. 

[pkVal pkInd] = max(interpWave(100:end));
pkInd = pkInd + 100 - 1;
%now we need to find the minimum following the peak

[troughVal troughInd] = min(interpWave(pkInd:end));
troughInd = troughInd + pkInd - 1;

peakTrough = (troughInd - pkInd)/300000;

if peakTrough < thresholdPeakValley
    funcOut = 1;
else
    funcOut = 0;
end

end
