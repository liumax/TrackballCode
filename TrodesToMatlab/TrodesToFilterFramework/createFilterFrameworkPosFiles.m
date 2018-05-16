function createFilterFrameworkPosFiles(dest,animalID,sessionNum)
%createFilterFrameworkPosFiles(dest,animalID,sessionNum)
%
%This function extracts position tracking information and saves data in the FilterFramework format. 
%It is assumed that there is at least one *.videoPositionTracking file in the current directory. 

%
%The function also tries to separate the data by epochs.  This 
%information is acquired from any .trodesComments files found in the 
%current directory, or by lookinbg for gaps in recording data if there is a '*.time' folder in the current directory
%conaining time information about the session (from
%extractTimeBinaryFile.m).
%
%
%dest -- the directory where the processed files should be saved for the
%animal
%animalID -- a string identifying the animal's id (appended to the
%beginning of the files).
%sessionNum -- the session number (in chronological order for the animal)


tFiles = dir(['*.videoPositionTracking']);
currDir = pwd;
sessionString = getTwoDigitNumber(sessionNum);

%Get the epoch boundaries.  This assumes that there is at least one
%*.trodesComments file in the current directory with "epoch start" and "epoch end" indicators 
epochList = getEpochs(1);  %assumes that there is at least a 1-second gap in data between epochs if no .trodesComments file is found
%epochs = readTrodesTaskFile();


posMatrix = [];
%There may be more than one file in the current directory.  If so, process
%them all, then sort by time.
for fileInd = 1:length(tFiles)
    disp(['Reading file: ',tFiles(fileInd).name]);
    offset = 0;
    tmpFileName = tFiles(fileInd).name;
    dotLoc = strfind(tmpFileName,'.');
    baseName = tmpFileName(1:dotLoc-1);
    tmpPosData = readTrodesExtractedDataFile(tmpFileName);
    if isfield(tmpPosData,'clockrate')
        clockrate = tmpPosData.clockrate;
    else
        clockrate = 30000;
    end

     %If a comments file exists, look to see if an offset was declared
    commentFile = dir([baseName,'.trodesComments']);
    if (~isempty(commentFile))
       fid = fopen(commentFile(1).name,'r');
    
        if (fid ~= -1)
        	
            while 1
                tline = fgetl(fid);
                if ~ischar(tline)
                    break
                end
               
                [fieldName, remainder] = strtok(tline);
                if isequal(fieldName,'offset')
                    remainder = remainder(2:end);
                    if ~isempty(str2num(remainder))                       
                        offset = str2num(remainder); 
                        disp(['Using offset: ', num2str(offset)]);
                    end
                end
                
            end
        end
    end
    
    %Create a matrix to hold the position data from the file
    tmpPosMatrix = zeros(length(tmpPosData.fields(1).data),5);
    tmpPosMatrix(:,1) = double(tmpPosData.fields(1).data + offset)/clockrate;
    for f = 2:length(tmpPosData.fields)
        tmpPosMatrix(:,f) = double(tmpPosData.fields(f).data);
    end
    posMatrix = [posMatrix; tmpPosMatrix];
    
end

%Sort the values by the timestamps.
posMatrix = sortrows(posMatrix,1);

fieldNames = {'time','x1','y1','x2','y2'};
rawpos{sessionNum} = [];
%check if epochs are defined, if so separate by epoch
if (~isempty(epochList))
    for e = 1:size(epochList,1)
        rawpos{sessionNum} = [];
        rawpos{sessionNum}{e} = [];
        epochPosMatrix = posMatrix(find((posMatrix(:,1) >= epochList(e,1)) & (posMatrix(:,1) < epochList(e,2))),:);
        
        for i=1:5
            rawpos{sessionNum}{e} = setfield(rawpos{sessionNum}{e},fieldNames{i},epochPosMatrix(:,i));
        end
        epochString = getTwoDigitNumber(e);
        cd(dest)
        save([animalID,'rawpos',sessionString,'-',epochString], 'rawpos');
        cd(currDir);
    end    
else
    for i=1:5
        rawpos{sessionNum} = setfield(rawpos{sessionNum},fieldNames{i},posMatrix(:,i));
    end
    cd(dest)
    save([animalID,'rawpos',sessionString], 'rawpos');
    cd(currDir);
end


%---------------------------------------------------

function numString = getTwoDigitNumber(input)
    
if (input < 10)
    numString = ['0',num2str(input)];
else
    numString = num2str(input);
end


    
