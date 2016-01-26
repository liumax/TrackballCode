%This is code that should group the use of mattias's functions for analysis
%of Trodes data.

%This extracts the binary file, and produces a folder with spike files of
%the format .dat. This folder will have one file per nTrode, which in my
%case, means 16 nTrodes.

%You must be in the directory with the actual recordings for this code to
%work.

fileName = '151217_ML151211A_R17_2400_noSound';

extractSpikeBinaryFiles(fileName);

%this should create matclust file folder
createAllMatclustFiles

%%
%This will then automatically generate the correct folder name and move to
%that folder.

folderName = strcat(fileName,'.spikes');
cd(folderName);


%%Once in this folder, we need to find all the files and process them one by
%one.

%this finds all contents of the folder, converts it to a table, and uses
%the table to generate cell arrays (direct conversion is messy).
spikeContents = dir;
spikeContents = struct2table(spikeContents);
dirContents = table2cell(spikeContents(:,4));
spikeNames = table2cell(spikeContents(:,1));

%this ugly for loop is for determining which dir components are folders,
%which we do not want. 
dirHolder = zeros(height(spikeContents),1);

for i = 1: height(spikeContents)
    dirHolder(i) = dirContents{i};
end
%finds directories in dir output, removes from spikeNames
dirMask = find(dirHolder==1);
spikeNames(dirMask) = [];
fieldNames = spikeNames;

for i = 1:length(fieldNames)
    y=strfind(fieldNames{i},'.');
    z= fieldNames{i}(y(1)+1:y(2)-1);
    fieldNames{i} = z;
end


%%
%this code will then use the readTrodesExtractedDataFile.m code to extract
%spike data from individual DAT files. 

cellData = struct;

for i = 1:length(spikeNames)
    cellData.(fieldNames{i}) = readTrodesExtractedDataFile(spikeNames{i});
end



