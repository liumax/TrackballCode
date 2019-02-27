




fnames = dir('Trial**');
numFiles = size(fnames,1);
fnames = {fnames.name};
fnames = fnames';

desigVal = zeros(numFiles,1);
for i = 1:numFiles
    tempFile = fnames{i};
    findSpace = strfind(tempFile,' ');
    findPer = strfind(tempFile,'.');
    desigVal(i) = str2num(tempFile(findSpace(end)+1:findPer-1));
    
end

for i = 1:numFiles
    v = VideoReader(fnames{i});
    vidHolder = zeros(352-136+1,996 - 578 + 1,v.NumberOfFrames,'int8');
%     vidHolder = uint8(vidHolder);
    for j = 1:v.NumberOfFrames
        tester = read(v,j);
        tester = mean(tester,3);
        tester = tester(136:352,578:996);
        tester = uint8(tester);
        vidHolder(:,:,j) = tester;
    end
    tester = read(v);
    v = [];
    tester = tester(136:352,578:996,:,:);
end

v = VideoReader('Trial     2.mpg')

tester = read(v);
tester = tester(136:352,578:996,:,:);