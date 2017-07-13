function [peakInfo, riseInfo, troughInfo] = findPhotoPeaks(t,corrData,thresh)
% find Peaks in photoSignal: actually finds peak magnitudes by looking at
% trough to peak distances. Does fine but could be improved for pulling out rise times. 

% peaks Amplitude is actual amplitude in dFoF
% rise amplitude is difference from trough to Peak (might change this to
        % rise location to peak??
% trough amplitude is actual amplitude at the trough


                             
mnPkDst = 2; % start of peak cannot be more than 2 sec before peak

% smooth signal 20 samples:
smSig = smooth(corrData,20);
dsSig = downsample(smSig,20); % downSampled
t_ds = downsample(t,20);


% IDEAL FILTER:
testSampleTS = timeseries(dsSig,t_ds);
mSig = mean(testSampleTS);
interval = [0.05 5];
idealF = idealfilter(testSampleTS,interval,'pass');
filtTime = idealF.Time;
filtSig = idealF.Data;
[locs] = findpeaks(double(filtSig));
pks = double(filtSig(locs.loc));
locs = locs.loc;


% find minimums:
sigInv = 1.01*max(filtSig) - filtSig;
[locsMn] = findpeaks(double(sigInv));
pksMn = double(filtSig(locsMn.loc));
locsMn = locsMn.loc;

% % Find derivative of filtered signal and only positive derivative:
% filtDerivSig = [0; diff(filtSig)];
% filtDerivSigPos = filtDerivSig; filtDerivSigPos(find(filtDerivSig<0)) = 0; % remove negatives
% filtDerivSig = interp1(filtTime,filtDerivSig,t); % interpolate to original size of corrData
% filtDerivSigPos = interp1(filtTime,filtDerivSigPos,t); % interpolate to original size of corrData



% THIS WILL BE THE ALGORITHM TO CHOOSE THROUGH PEAKS:
count = 1;
rmPkCnt = 0;
for i = 1:length(pks)
    mnIdx  = find(locsMn < locs(i), 1, 'last'); % find index of last trough

    % if no trough exists before peak, skip to next peak
    if ~any(mnIdx)
        continue;
    end
    
    if filtTime(locsMn(mnIdx)) + mnPkDst < filtTime(locs(i))
        rmPkCnt = rmPkCnt+1;
%         disp(['Peak ', num2str(i), ...
%             ' excluded because trough-peak dist too great: ', ...
%             num2str(filtTime(locs(i)) - filtTime(locsMn(mnIdx)))]);
        [mnVal,mnIdx] = find(filtTime > filtTime(locs(i)) - mnPkDst,1,'first');
    end
    
    % Find the max slope between the trough and peak
    [pkSlope,pkSlopeInd] = max(diff(filtSig(locsMn(mnIdx):locs(i))));
    rIdx = locsMn(mnIdx) + pkSlopeInd; % peak slope index

    
    [pkRise,pkRiseInd] = max(diff(diff(filtSig(locsMn(mnIdx):rIdx)))); % 2nd derivative (trough to rise)
    r2Idx = locsMn(mnIdx) + pkRiseInd;
    
    if ~any(r2Idx)
        continue;
    end
    
    % Store only peaks that exceed a threshold
    if (pks(i)-filtSig(locsMn(mnIdx))) > thresh
        
        peakInfo.t(count) = filtTime(locs(i));
        peakInfo.sampleNum(count) = find(t >= filtTime(locs(i)), 1, 'first');
        peakInfo.amp(count) = corrData(peakInfo.sampleNum(count)); % Actual amplitude at peak

        
        
        troughInfo.t(count) = filtTime(locsMn(mnIdx));
        troughInfo.sampleNum(count) = find(t >= filtTime(locsMn(mnIdx)), 1, 'first');
        troughInfo.amp(count) = corrData(troughInfo.sampleNum(count)); % actual amplitude at trough
        
        riseInfo.t(count) = filtTime(r2Idx);
        riseInfo.sampleNum(count) = find(t >= filtTime(r2Idx), 1, 'first');
        riseInfo.amp(count) = peakInfo.amp(count) - troughInfo.amp(count); % RELATIVE AMPITUDES (Peak-Trough)

        
        
        count = count+1;
    end
end    





% Remove outliers: (get rid of anything >4 std below minimum)
zthresh = -4;
[z] = zscore(corrData(troughInfo.sampleNum))';
rmInds = find(z<zthresh | riseInfo.amp<0);
peakInfo.t(rmInds) = [];
peakInfo.sampleNum(rmInds) = [];
peakInfo.amp(rmInds) = [];
riseInfo.t(rmInds) = [];
riseInfo.sampleNum(rmInds) = [];
riseInfo.amp(rmInds) = [];
troughInfo.t(rmInds) = [];
troughInfo.sampleNum(rmInds) = [];
troughInfo.amp(rmInds) = [];


% Throw warning if too many peaks were removed:
disp([num2str(100*rmPkCnt/length(peakInfo.t),'%1.1f'),'% PEAKS REMOVED',...
    ': trough-peak greater than ',num2str(mnPkDst),'s'])
if (100*rmPkCnt/length(peakInfo.t))>5
    warning('TOO MANY PEAKS REMOVED: CHECK FILE!!!')
end
        
    
% Uncomment to see peaks and amplitudes:
figure
subplot(2,1,1)
hold on
plot(t,corrData,'k')
hold on
axis([peakInfo.t(1) peakInfo.t(end) -.02 0.3])
plot(peakInfo.t,corrData(peakInfo.sampleNum),'b.')
plot(riseInfo.t,corrData(riseInfo.sampleNum),'g.')
plot(troughInfo.t,corrData(troughInfo.sampleNum),'r.')
subplot(2,1,2)
hold on
hist(riseInfo.amp,200)
xlabel('Amplitudes')
ylabel('Frequency')
title('Trough to peak amplitude distribution')




