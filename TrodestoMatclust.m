
%this is the filename for the .rec file. Please put it in the same dir as
%the config file if you decide to do post-hoc referencing etc.
fileName = 'pairingTester';

%this extracts timestamps
extractTimeBinaryFile(fileName)

%this extracts spikes with respect to any changes in the config file
extractSpikeBinaryFiles(fileName)

%this generates matclust files.
createAllMatclustFiles

%this will extract LFP data
extractLFPBinaryFiles(fileName, 0);

%this will extract DIO data (inputs and outputs)
extractDioBinaryFiles(fileName);

%finds all files related to the recording
dirStruct = dir; %generates directory
dirNames = {dirStruct.name}; %pulls file/folder names
allMatches = strfind(dirNames,fileName); %finds things matching target file name
allMatches = find(~cellfun(@isempty,allMatches)==1); %converts empty cells to logical and finds where it is true
dirConf = find(arrayfun(@(y)all(y.isdir==1),dirStruct)==1); %finds all things that are directories
moveDir = intersect(allMatches,dirConf);

%generates directory to put everything into!
mkdir(fileName(1:end-2)); %makes directory
for i = 1:size(moveDir,1) %moves all folders in!!
    movefile(dirNames{moveDir(i)},fileName(1:end-2));
end

currDir = pwd; %finds current path
newDir = strcat(currDir,'\',fileName(1:end-2)); %generates path for new directory
cd(newDir)
%finds directory, identifies matclust folder!
dirStruct = dir;
dirNames = {dirStruct.name};
allMatches = strfind(dirNames,'matclust');
allMatches = find(~cellfun(@isempty,allMatches)>0);
%finds path for matclust folder
currDir = pwd;
newDir = strcat(currDir,'\',dirNames{allMatches});
%changes directory to matclust directory
cd(newDir)
%open matclust
% matclust


% %this will extract continuous data. as of 160216 doesnt work. 
% % !"C:\Trodes\exportLFP" -rec  160129_ML150108A_L17_3000_toneFinder.rec -userefs 0
% trodesPath = which('trodes_path_placeholder.m');
% trodesPath = fileparts(trodesPath);
% 
% contDataCommand = ['!"',fullfile(trodesPath,'exportLFP'),'"', ' -rec ',fileName,'.rec',...
%     ' -usespikefilters 0 -lowpass -1 -highpass -1 -userefs 0 -outputrate 30000'];
% 
% eval(contDataCommand)