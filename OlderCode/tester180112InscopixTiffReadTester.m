





fileName = '180105_ML171117A_GMC.tif'









tester = imfinfo(fileName);
numFrames = length(tester);


for i = 1:length(alignTimes)
    i
    %get pre-frames
    for j = 1:10
        tiff_stack = imread('180105_ML171117A_GMC.tif', [alignTimes(i) - j]) ; 
        stackStore(:,:,j) = tiff_stack;
    end
    meanPre = mean(stackStore,3);
    for j = 1:10
        tiff_stack = imread('180105_ML171117A_GMC.tif', [alignTimes(i) + j]) ; 
        stackStore(:,:,j) = tiff_stack;
    end
    meanPost = mean(stackStore,3);
    fullStore(:,:,i) = meanPost-meanPre;
    fullRawStore{i,1} = meanPre;
    fullRawStore{i,2} = stackStore;
end


photoAverages = zeros(size(fullStore,1),size(fullStore,2),numFreqs,numDBs);
photoRasterStore = cell(numFreqs,numDBs);
for ind = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(ind),3));
    for k = 1:numDBs
        %find all trials of the particular setting
        targetFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == subUniqueDB(k));
        %pull and average these traces
        tempHolder = fullStore(:,:,targetFinder);
        photoRasterStore{ind,k} = tempHolder;
        tempHolder = mean(tempHolder,3);
        photoAverages(:,:,ind,k) = tempHolder;
    end
end

fullMin = min(min(min(photoAverages)));
fullMax = max(max(max(photoAverages)));
for i = 1:14
figure
imagesc(photoAverages(:,:,i),[fullMin,fullMax])
end
