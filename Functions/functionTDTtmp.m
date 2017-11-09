%This function is meant to perform extraction and then save as a tmp file,
%in case further analysis code fails. 

function [filtSig1,filtSig2,traceDF,traceTiming,t_ds,newSmoothDS,targetPeaks,data] = functionTDTtmp(fileName,zToggle);

%first, look for tmp file! first, we want to pull all files!
folderFiles = what;
folderFiles= folderFiles.mat;

tmpName = strcat(fileName,'TDTTMP.mat');
[findString] = functionCellStringFind(folderFiles,tmpName);
disp('LOOKING FOR TDT TMP FILE')
if findString %if there is a tmp file!
    disp('TDT TMP FILE FOUND! LOADING')
    load(folderFiles{findString})
else
    disp('NO TMP FOR TDT DATA, EXTRACTING...')
    %load file
    data = load(strcat(fileName,'.mat'));
    data=data.data;

    [filtSig1,filtSig2,traceDF,traceTiming] = functionPhotometryRawExtraction(data);

    %pull peaks 170616 This appears to have problem: built for 2016 matlab, has
    %additional functionality for peak finding.
    try
        if zToggle == 0
            [t_ds,newSmoothDS,targetPeaks] = functionPhotoPeakProcess(traceTiming,filtSig1,0.005);
        else
            [t_ds,newSmoothDS,targetPeaks] = functionPhotoPeakProcessZ(traceTiming,(filtSig1-mean(filtSig1))/std(filtSig1),0.1);
        end
    catch
        error('Peak Detection Failed')
%         targetPeaks = [];
%         newSmoothDS = [];
%         t_ds = [];
    end
    tmpName = strcat(fileName,'TDTTMP.mat');
    save(tmpName,'filtSig1','filtSig2','traceDF','traceTiming','t_ds','newSmoothDS','targetPeaks','data');
    disp('TDT DATA SAVED AS TMP')
end

end