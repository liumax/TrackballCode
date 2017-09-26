%This code is meant to go through tuning data and try to pick out reliable
%responses that can then be used for behavior!

%load target files


testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);

[findString] = functionCellStringFind(testNames,'Sound');
testNames(findString) = [];

[findString] = functionCellStringFind(testNames,'Tuning');
testNames = testNames(findString);


numFiles = length(testNames);

for bigInd = 1:numFiles
    
    load(testNames{bigInd})
    
    %lets just plot things?
    figure
    plot(squeeze(s.Processed.PhotoAverages(:,:,3)))
    
    %now lets go through the data and extract the shits. only looking at
    %maximum amplitude for now
    numFreqs = size(s.Processed.PhotoStore,1);
    
    for i = 1:numFreqs
        tempData = s.Processed.PhotoStore{i,3};
        baseVal = min(tempData(140:150,:));
        peakVal = max(tempData(150:170,:));
        peakMag = peakVal - baseVal;
        magStore(:,i) = peakMag;
        meanStore(i) = mean(peakMag);
        stdStore(i) = std(peakMag);
        corrStore{i} = tempData - repmat(baseVal,382,1);
    end
    
    figure
    plot(meanStore,'LineWidth',2)
    hold on
    plot(meanStore - stdStore)
    plot(meanStore + stdStore)
end




