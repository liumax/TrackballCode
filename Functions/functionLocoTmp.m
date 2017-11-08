%This function is meant to perform extraction and then save as a tmp file,
%in case further analysis code fails. 

function [locoData] = functionLocoTmp(fileName,portStates,locoTimeStep);

%first, look for tmp file! first, we want to pull all files!
folderFiles = what;
folderFiles= folderFiles.mat;

tmpName = strcat(fileName,'LocoTMP.mat');

[findString] = functionCellStringFind(folderFiles,tmpName);
disp('LOOKING FOR LOCO TMP FILE')

if findString %if there is a tmp file!
    disp('LOCO TMP FILE FOUND! LOADING')
    load(folderFiles{findString})
else
    disp('NO TMP FOR LOCO DATA, EXTRACTING...')
    [locoData] = functionMBEDrotary(portStates.inStates(:,4),portStates.inStates(:,5),portStates.tStamps/1000,locoTimeStep);
    save(tmpName,'locoData');
    disp('LOCODATA SAVED AS TMP')
end

end