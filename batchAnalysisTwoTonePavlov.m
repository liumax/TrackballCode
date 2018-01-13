%This will be a wrapper function for the photometry analysis code. 

testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames(findString) = [];
[findString] = functionCellStringFind(testNames,'Pavlov');
testNames = testNames(findString);
[findString] = functionCellStringFind(testNames,'TMP');
testNames(findString) = [];

numFiles = length(testNames);
badInd = 1
for bigInd = 1:numFiles
    fileName = testNames{bigInd}(1:end-4);
    try
        [s] = analysisTwoTonePavlovPhotometry(fileName);
        s = [];
    catch
        badStore(badInd) = bigInd;
        badInd = badInd + 1;
        disp('CODE FAILURE, SAVING FILE')
        diary off
    end
end