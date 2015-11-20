dirName = 'E:\tempAcuteRecordings\Processed\ML151104A';
files = dir(fullfile(dirName,'*.nex'));
files = {files.name};

batchDataStructure = struct;

for i=1:length(files)
    fileName = files{i};
    [nexFile] = readNexFile(fileName);
    [s] = auditoryStatsAnalysis(nexFile,fileName);
    x = fileName;
    y = strfind(x,'.');
    x = x(1:y-1);
    batchDataStructure.(x) = s;
end