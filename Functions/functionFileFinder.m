function [targetFiles] = functionFileFinder(subFoldersCell,folderTarget,fileTarget);
%This code is meant to find a folder (folderTarget) from a range of folders 
%(subFoldersCell), and search that folder for specific contents 
%(fileTarget), and make a cell array of said files.
folderFinder = strfind(subFoldersCell,folderTarget);%finds DIO folder
folderFinder = find(~cellfun(@isempty,folderFinder)); %determines empty cells, finds index for target cell
targetFolderName = subFoldersCell{folderFinder}; %pulls out the name for the DIO folder
targetFolderSearch = dir(targetFolderName);%pulls dir from DIO folder
targetFileNames = {targetFolderSearch.name}';%pulls names section
targetFileFinder = strfind(targetFileNames,fileTarget); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name