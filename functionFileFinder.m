function [targetFiles] = functionFileFinder(subFoldersCell,folderTarget,fileTarget);
%find DIO folder and D1 file for analysis
folderFinder = strfind(subFoldersCell,folderTarget);%finds DIO folder
folderFinder = find(~cellfun(@isempty,folderFinder)); %determines empty cells, finds index for target cell
targetFolderName = subFoldersCell{folderFinder}; %pulls out the name for the DIO folder
targetFolderSearch = dir(targetFolderName);%pulls dir from DIO folder
targetFileNames = {targetFolderSearch.name}';%pulls names section
targetFileFinder = strfind(targetFileNames,fileTarget); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name