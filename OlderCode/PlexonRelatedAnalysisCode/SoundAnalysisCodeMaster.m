dirName = 'E:\626CDatStimAnalysisForF30';
% dirName = 'E:\tempAcuteRecordings\Processed\ML151104A';
files = dir(fullfile(dirName,'*.nex'));
files = {files.name};

batchDataStructure = struct;

for i=1:length(files)
    fileName = files{i};
    [nexFile] = readNexFile(fileName);
    %note that you need to be in the directory with the files for
    %readNexFile to work!
    [s] = auditoryStatsAnalysis(nexFile,fileName);
    x = fileName;
    x = genvarname(x);
    batchDataStructure.(x) = s;
%     batchDataStructure.(x) = 0;
    [output] = auditoryGraphing(nexFile,fileName);
    batchDataStructure.(x).Plotting = output;
    [statsGraph] = auditoryStatGraph(batchDataStructure,x);
    batchDataStructure.(x).ResponsesStats = statsGraph;
end