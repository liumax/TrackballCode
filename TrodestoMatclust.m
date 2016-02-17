
%this is the filename for the .rec file. Please put it in the same dir as
%the config file if you decide to do post-hoc referencing etc.
fileName = '160129_ML150108A_L17_3000_toneFinder';

%this extracts timestamps
extractTimeBinaryFile(fileName)

%this extracts spikes with respect to any changes in the config file
extractSpikeBinaryFiles(fileName)

%this generates matclust files.
createAllMatclustFiles

%this will extract LFP data
extractLFPBinaryFiles(fileName, 0);

% %this will extract continuous data. as of 160216 doesnt work. 
% !"C:\Trodes\exportLFP" -rec  160129_ML150108A_L17_3000_toneFinder.rec -userefs 0
% contDataCommand = ['!"',fullfile(trodesPath,'exportLFP'),'"', ' -rec ',fileName,'.rec',...
%     ' -usespikefilter 0 -lowpass -1 -highpass -1 -userefs 0 -outputrate 30000'];