%180629 Testing code to look at video frames
%goal here is to extract video frames at peaks of "oscillating" activity.

Wn = [0.5*2/100 5*2/100]; % %confirmed 180601 that this is correct
%     formula. do 2 x desired Hz /sampling rate
n = 1000; % 1000th order filter (slower? but 100-order was too low)
b = fir1(n, Wn); 
filtData = filtfilt(b,1,s.nt2_56.FineSession);
figure
plot(filtData)
findRun = find(s.RotaryData.NewSimp == 1);
fundRunTimes = (s.RotaryData.Velocity(findRun,1));
%find peaks
peakLocs = findpeaks(filtData);
peakLocs = peakLocs.loc;

peakMags = filtData(peakLocs);
%find peak mags above 0.2
bigPeaks = find(peakMags > 0.2);
bigPeakTimes = s.nt2_56.FineSessionVector(peakLocs(bigPeaks));
%find matching frames?
matchTimes = [];
for i = 1:length(bigPeakTimes)
    matchTimes(i) = find(s.CameraTimes < bigPeakTimes(i),1,'last');
end
