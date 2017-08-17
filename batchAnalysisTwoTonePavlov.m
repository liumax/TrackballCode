%This will be a wrapper function for the photometry analysis code. 

testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames(findString) = [];

numFiles = length(testNames);

for bigInd = 1:numFiles
    fileName = testNames{bigInd}(1:end-4);
    [s] = analysisTwoTonePavlovPhotometry(fileName);
    s = [];
end