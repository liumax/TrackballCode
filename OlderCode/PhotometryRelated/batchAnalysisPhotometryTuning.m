


%This will be a wrapper function for the photometry analysis code. 

testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames(findString) = [];

[findString] = functionCellStringFind(testNames,'Sound');
testNames(findString) = [];

[findString] = functionCellStringFind(testNames,'TMP');
testNames(findString) = [];

% [findString] = functionCellStringFind(testNames,'Tuning');
% testNames = testNames(findString);

numFiles = length(testNames);

for bigInd = 1:numFiles
    fileName = testNames{bigInd}(1:end-4);
    [s] = analysisPhotometryTuning(fileName);
    
end