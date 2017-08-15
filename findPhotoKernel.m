function P = findPhotoKernel(t,photoSig,evTimes)
% Based on Parker/Witten (2016) method for extracting 
% t is time vector in ms
% photoSig is photoSignal
% evTimes is a cell array indicating the timing of events that the user
% wants. You can use any number of events to put into the model

% 
% NOTE: this has a hard-coded downsamping factor of 50, which corresponds
% roughly to 50ms sampling for the kernel
% CHD: 17.07.07 MHL170814 EDITED: Now automatically detects, and
% downsamples based on detected frequency of sampling. Rounds down
% down-sample factor in favor of more points rather than less.

% Downsample signal:
sampleRate = mean(diff(t));
sampleRatio = floor(0.05/sampleRate);
d = downsample([t' photoSig],sampleRatio); % want to target downsampling to 50 ms bins

% d = [t' photoSig];

% Create time vector for Kernel:
dt = mean(diff(d(:,1)));
a = fliplr([0:-dt:-1]);
b = dt:dt:2;
tVec = [a b];

% Create event Vectors:
for kerNum = 1:size(evTimes,2)
    evVec{kerNum} = zeros(size(d,1),1); % initialize
    for evNum = 1:size(evTimes{kerNum},1)
        idx = find(abs(d(:,1)-evTimes{kerNum}(evNum))== ...
        min(abs(d(:,1)-evTimes{kerNum}(evNum))));
        evVec{kerNum}(idx) = 1;
    end
end
    
% Create Predictor Matrix:
X = [];
for kerNum = 1:size(evVec,2)
    for i = 1:length(tVec)
        evConv{kerNum}(:,i) = circshift(evVec{kerNum},i-length(a));
    end
    X = [X evConv{kerNum}];
end

% FIT MODEL:
[b dev stats] = glmfit(X,zscore(d(:,2)));


% Store data:
P.t = tVec;
P.const = b(1); % constant term 
stIdx = 2; % 
for kerNum = 1:size(evVec,2)
    edIdx = stIdx+length(tVec)-1;
    P.kernel{kerNum}.b = b(stIdx:edIdx);
    P.kernel{kerNum}.pVal = stats.p(stIdx:edIdx);
    stIdx = edIdx+1;
end


