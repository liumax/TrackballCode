%This function is meant to perform extraction and then save as a tmp file,
%in case further analysis code fails. 

function [trialStates, portStates, trialParams,inputPhot,outputPhot,photoToggle] = functionMBEDtmp(fileName);

%first, look for tmp file! first, we want to pull all files!
folderFiles = what;
folderFiles= folderFiles.mat;
%set TMP name
tmpName = strcat(fileName,'MBEDTMP.mat');

[findString] = functionCellStringFind(folderFiles,tmpName);
disp('LOOKING FOR MBED TMP FILE')
if findString %if there is a tmp file!
    disp('MBED TMP FILE FOUND! LOADING')
    load(folderFiles{findString})
else
    disp('NO MBED TMP FILE, EXTRACTING...')
    
    [trialStates, portStates, trialParams] = maxTrialVariablesLickingTask(fileName);

    %we want to look at instates. Port 1 is TDT for photometry, port 2 is
    %NOLDUS

    inputPhot = [portStates.tStamps',portStates.inStates(:,8)];
    outputPhot = [portStates.tStamps',portStates.outStates(:,8)];

    %eliminate duplicate values!
    try
        [inputPhot] = functionSignalDuplicateElim(inputPhot,2);
        photoToggle = 0;
        disp('Duplicate Elimination for Inputs Successful!')
    catch
        disp('NO INPUTS, SWITCH TO USING PHOTOMETRY')
        photoToggle = 1;
    end
    [outputPhot] = functionSignalDuplicateElim(outputPhot,2);
    %save temporary file so that i can save time later on!
    
    save(tmpName,'trialStates','portStates','trialParams','inputPhot','outputPhot','photoToggle');
    disp('SAVED TMP FILE FOR MBED')
end

end