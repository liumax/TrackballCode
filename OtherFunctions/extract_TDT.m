% Script to extract TDT files and organize


startPath = 'C:\TDT\Synapse\Tanks';
[dirName] = uigetdir(startPath);
data = TDT2mat(dirName);
[fileName,pathName] = uiputfile;
save([pathName,fileName],'data')
